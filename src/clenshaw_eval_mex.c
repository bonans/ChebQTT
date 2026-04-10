/*
 * CLENSHAW_EVAL_MEX  Evaluate sum of Chebyshev polynomials via Clenshaw.
 *
 *   Y = CLENSHAW_EVAL_MEX(X, V, A, B)
 *
 *   X : (N x 1) evaluation points
 *   V : (m+1 x M) coefficient matrix (ascending degree k=0,1,...,m)
 *   A, B : scalars defining the domain [A, B]
 *
 *   Y : (N x M) output,  Y(i,j) = sum_{k=0}^{m} V(k+1,j) * T_k( x_mapped(i) )
 *
 *   Uses Clenshaw recurrence so only O(N*M) memory.
 *
 *   Optimizations (compared to v1):
 *     - Two-step unrolled Clenshaw (halves loop iterations, matches Chebfun)
 *     - Cache-blocked evaluation (bk1, bk2, t2 fit in L1 cache)
 *     - Precomputed t2 per block (eliminates per-step recomputation)
 *     - restrict pointers for alias-free optimization
 *     - Signed loop variables for better compiler vectorization
 *
 *   Compile:
 *     mex COPTIMFLAGS="-O3 -march=native -ffast-math -DNDEBUG" clenshaw_eval_mex.c
 */

#include "mex.h"
#include <string.h>

/* Block size: working set = bk1(BS) + bk2(BS) + t2(BS) = 3*BS*8 bytes.
 * BS=1024 → 24KB, fits comfortably in L1 cache (typically 32-48KB). */
#define BS 1024

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --- Input validation --- */
    if (nrhs != 4)
        mexErrMsgIdAndTxt("clenshaw_eval_mex:nrhs",
                          "Four inputs required: X, V, A, B.");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("clenshaw_eval_mex:nlhs",
                          "At most one output.");

    /* X: N x 1 */
    const double * restrict X = mxGetPr(prhs[0]);
    const mwSize N = mxGetNumberOfElements(prhs[0]);

    /* V: (m+1) x M */
    const double * restrict V = mxGetPr(prhs[1]);
    const mwSize L = mxGetM(prhs[1]);   /* m+1 */
    const mwSize M = mxGetN(prhs[1]);

    /* Domain */
    const double a = mxGetScalar(prhs[2]);
    const double b = mxGetScalar(prhs[3]);
    const double scale = 2.0 / (b - a);
    const double shift = -(a + b) / (b - a);

    /* Output: N x M */
    plhs[0] = mxCreateDoubleMatrix(N, M, mxREAL);
    double * restrict Y = mxGetPr(plhs[0]);

    const long m = (long)L - 1;   /* signed for clean loops */
    mwSize j, ib;

    if (m == 0) {
        /* Constant polynomial: Y(i,j) = V(1,j) */
        for (j = 0; j < M; j++) {
            const double c0 = V[j * L];
            double * restrict Yj = Y + j * N;
            for (mwSize i = 0; i < N; i++)
                Yj[i] = c0;
        }
        return;
    }

    /* --- Two-step Clenshaw recurrence with cache blocking ---
     *
     * Standard Clenshaw for T_k on [-1,1]:
     *   b_{m+1} = 0,  b_{m+2} = 0
     *   b_k = c_k + 2t * b_{k+1} - b_{k+2}   (k = m, m-1, ..., 1)
     *   y   = c_0 + t * b_1 - b_2
     *
     * Two-step unrolling (Chebfun style) processes two steps per iteration:
     *   bk2 = c_k   + 2t * bk1 - bk2    → b_k
     *   bk1 = c_{k-1} + 2t * bk2 - bk1  → b_{k-1}
     *
     * Block over points to keep bk1[], bk2[], t2[] in L1 cache.
     */

    /* Stack-allocated workspace for one block */
    double bk1[BS], bk2[BS], t2[BS];

    for (j = 0; j < M; j++) {
        const double * restrict Vj = V + j * L;  /* column j of V */

        for (ib = 0; ib < N; ib += BS) {
            const mwSize bs = (ib + BS <= N) ? BS : (N - ib);
            const double * restrict Xb = X + ib;
            double * restrict Yb = Y + j * N + ib;
            mwSize i;
            long k;

            /* Precompute t2 = 2*mapped_x for this block */
            for (i = 0; i < bs; i++) {
                t2[i] = 2.0 * (Xb[i] * scale + shift);
                bk1[i] = 0.0;
                bk2[i] = 0.0;
            }

            /* Two-step Clenshaw: from k = m down to 2 (step -2) */
            for (k = m; k >= 2; k -= 2) {
                const double ck  = Vj[k];
                const double ck1 = Vj[k - 1];
                for (i = 0; i < bs; i++) {
                    /* b_k = c_k + 2t * b_{k+1} - b_{k+2} */
                    double new_bk2 = ck + t2[i] * bk1[i] - bk2[i];
                    /* b_{k-1} = c_{k-1} + 2t * b_k - b_{k+1} */
                    bk1[i] = ck1 + t2[i] * new_bk2 - bk1[i];
                    bk2[i] = new_bk2;
                }
            }

            /* Handle remaining step for odd degree */
            if (m % 2 == 1) {
                const double c1 = Vj[1];
                for (i = 0; i < bs; i++) {
                    double tmp = bk1[i];
                    bk1[i] = c1 + t2[i] * bk1[i] - bk2[i];
                    bk2[i] = tmp;
                }
            }

            /* Final step: y = c_0 + t * b_1 - b_2
             * where t = mapped_x (NOT 2*t), so t = t2/2 */
            {
                const double c0 = Vj[0];
                for (i = 0; i < bs; i++)
                    Yb[i] = c0 + 0.5 * t2[i] * bk1[i] - bk2[i];
            }
        }
    }
}
