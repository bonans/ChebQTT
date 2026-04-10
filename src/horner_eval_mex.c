/*
 * HORNER_EVAL_MEX  Evaluate sum of monomials via Horner's scheme.
 *
 *   Y = HORNER_EVAL_MEX(X, V)
 *
 *   X : (N x 1) evaluation points
 *   V : (m+1 x M) coefficient matrix (ascending degree k=0,1,...,m)
 *
 *   Y : (N x M) output,  Y(i,j) = sum_{k=0}^{m} V(k+1,j) * X(i)^k
 *
 *   Uses Horner's method so only O(N*M) memory (no Vandermonde matrix).
 *
 *   Optimizations:
 *     - Cache-blocked evaluation (block fits in L1 cache)
 *     - restrict pointers for alias-free optimization
 *     - Signed loop variables for better compiler vectorization
 *     - Two-step loop unrolling to reduce loop overhead
 *
 *   Compile:
 *     mex COPTIMFLAGS="-O3 -march=native -ffast-math -DNDEBUG" horner_eval_mex.c
 */
#include "mex.h"
#include <string.h>

/* Block size: working set = Xb(BS) + Yb(BS) = 2*BS*8 bytes.
 * BS=2048 → 32KB, fits in L1 cache (typically 32-48KB). */
#define BS 2048

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --- Input validation --- */
    if (nrhs != 2)
        mexErrMsgIdAndTxt("horner_eval_mex:nrhs",
                          "Two inputs required: X, V.");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("horner_eval_mex:nlhs",
                          "At most one output.");

    /* X: N x 1 */
    const double * restrict X = mxGetPr(prhs[0]);
    const mwSize N = mxGetNumberOfElements(prhs[0]);

    /* V: (m+1) x M */
    const double * restrict V = mxGetPr(prhs[1]);
    const mwSize L = mxGetM(prhs[1]);   /* m+1 */
    const mwSize M = mxGetN(prhs[1]);
    const long m = (long)L - 1;          /* signed for clean loops */

    /* Output: N x M */
    plhs[0] = mxCreateDoubleMatrix(N, M, mxREAL);
    double * restrict Y = mxGetPr(plhs[0]);

    mwSize j, ib;

    for (j = 0; j < M; j++) {
        const double * restrict Vj = V + j * L;  /* column j of V */

        /* --- Blocked Horner evaluation ---
         * Process BS points at a time so Xb and Yb stay in L1 cache
         * across all m Horner steps. */
        for (ib = 0; ib < N; ib += BS) {
            const mwSize bs = (ib + BS <= N) ? BS : (N - ib);
            const double * restrict Xb = X + ib;
            double * restrict Yb = Y + j * N + ib;
            mwSize i;
            long k;

            /* Initialize Yb = c_m */
            const double cm = Vj[m];
            for (i = 0; i < bs; i++)
                Yb[i] = cm;

            /* Two-step unrolled Horner: process 2 degrees per iteration */
            for (k = m - 1; k >= 1; k -= 2) {
                const double ck  = Vj[k];
                const double ck1 = Vj[k - 1];
                for (i = 0; i < bs; i++) {
                    Yb[i] = ck + Xb[i] * Yb[i];
                    Yb[i] = ck1 + Xb[i] * Yb[i];
                }
            }
            /* Handle remaining step if m is odd (odd number of Horner steps) */
            if (m % 2 == 1) {
                const double c0 = Vj[0];
                for (i = 0; i < bs; i++)
                    Yb[i] = c0 + Xb[i] * Yb[i];
            }
        }
    }
}
