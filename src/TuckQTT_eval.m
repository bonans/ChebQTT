function vals = TuckQTT_eval(TuckQTT, x, method)
%TUCKQTT_EVAL Evaluate a function in QTT-Tucker format at given points.
%   VALS = TUCKQTT_EVAL(TUCKQTT, X) evaluates the QTT-Tucker tensor at
%   the points specified by X using nearest neighbor interpolation.
%
%   VALS = TUCKQTT_EVAL(TUCKQTT, X, METHOD) uses the specified interpolation
%   method: 'nearest' (default) or 'linear' for trilinear interpolation.
%
%   Input:
%     TUCKQTT - Struct in QTT-Tucker format with fields:
%       .core  - Core tensor in Tucker format (full format)
%       .cols  - Joint QTT for polynomials in first mode
%       .rows  - Joint QTT for polynomials in second mode
%       .tubes - Joint QTT for polynomials in third mode
%       .dom   - Domain of the function (3x2 matrix)
%       .d     - Number of quantics dimensions (3x1 vector)
%     X - Points at which to evaluate, tensor of size n1 x ... x nk x 3
%         Each slice X(i1,...,ik,:) is a point (x,y,z) in 3D space.
%     METHOD - 'nearest' or 'linear' (default: 'nearest')
%
%   Output:
%     VALS - Evaluated function values, tensor of size n1 x ... x nk
%
%   See also: CHEBTUCK2TUCKQTT, POLY2QTT

    if nargin < 3
        method = 'nearest';
    end

    dom = TuckQTT.dom;
    d = TuckQTT.d;
    core = TuckQTT.core;

    % Get Tucker core dimensions
    R = size(core);

    % Grid sizes
    n = 2 .^ d;
    % Grid step sizes
    h = (dom(:,2) - dom(:,1)) ./ n;

    % Reshape x to N x 3
    sz = size(x);
    N = prod(sz(1:end-1));
    x = reshape(x, N, 3);

    if strcmp(method, 'nearest')
        % Nearest neighbor interpolation
        % Find the nearest grid index for each point
        idx = zeros(N, 3);
        for dim = 1:3
            % Grid points are at: dom(dim,1) + (i - 0.5) * h(dim) for i = 1, ..., n(dim)
            % Find nearest: round((x - dom(dim,1)) / h(dim) + 0.5)
            idx(:, dim) = round((x(:, dim) - dom(dim,1)) / h(dim) + 0.5);
            % Clamp indices to be within [1, n(dim)]
            idx(:, dim) = max(1, min(n(dim), idx(:, dim)));
        end

        % Convert linear indices to binary multi-indices
        subs1 = lin2sub_all2(d(1), idx(:,1));  % N x d(1)
        subs2 = lin2sub_all2(d(2), idx(:,2));  % N x d(2)
        subs3 = lin2sub_all2(d(3), idx(:,3));  % N x d(3)

        % Batch-evaluate QTT tensors at all required indices
        F1 = qtt_batch_eval(TuckQTT.cols, subs1);   % N x r1
        F2 = qtt_batch_eval(TuckQTT.rows, subs2);   % N x r2
        F3 = qtt_batch_eval(TuckQTT.tubes, subs3);  % N x r3

        % Contract Tucker core with factors
        vals = contract_tucker(core, F1, F2, F3, N, R);

    elseif strcmp(method, 'linear')
        % Trilinear interpolation
        % Find the lower corner indices for each point
        idx_lo = zeros(N, 3);
        t = zeros(N, 3);  % Interpolation weights
        for dim = 1:3
            % Grid points are at: dom(dim,1) + (i - 0.5) * h(dim) for i = 1, ..., n(dim)
            % Continuous index: (x - dom(dim,1)) / h(dim) + 0.5
            cont_idx = (x(:, dim) - dom(dim,1)) / h(dim) + 0.5;
            idx_lo(:, dim) = floor(cont_idx);
            % Clamp lower index to [1, n(dim)-1] to ensure upper index is valid
            idx_lo(:, dim) = max(1, min(n(dim)-1, idx_lo(:, dim)));
            % Interpolation weight: distance from lower grid point (normalized by h)
            t(:, dim) = cont_idx - idx_lo(:, dim);
            % Clamp t to [0, 1] for boundary cases
            t(:, dim) = max(0, min(1, t(:, dim)));
        end
        idx_hi = idx_lo + 1;

        % Convert all 8 corner indices to binary multi-indices
        subs1_lo = lin2sub_all2(d(1), idx_lo(:,1));
        subs1_hi = lin2sub_all2(d(1), idx_hi(:,1));
        subs2_lo = lin2sub_all2(d(2), idx_lo(:,2));
        subs2_hi = lin2sub_all2(d(2), idx_hi(:,2));
        subs3_lo = lin2sub_all2(d(3), idx_lo(:,3));
        subs3_hi = lin2sub_all2(d(3), idx_hi(:,3));

        % Batch-evaluate QTT tensors at lower and upper indices
        F1_lo = qtt_batch_eval(TuckQTT.cols, subs1_lo);
        F1_hi = qtt_batch_eval(TuckQTT.cols, subs1_hi);
        F2_lo = qtt_batch_eval(TuckQTT.rows, subs2_lo);
        F2_hi = qtt_batch_eval(TuckQTT.rows, subs2_hi);
        F3_lo = qtt_batch_eval(TuckQTT.tubes, subs3_lo);
        F3_hi = qtt_batch_eval(TuckQTT.tubes, subs3_hi);

        % Interpolation weights
        t1 = t(:, 1); t2 = t(:, 2); t3 = t(:, 3);
        w1_lo = 1 - t1; w1_hi = t1;
        w2_lo = 1 - t2; w2_hi = t2;
        w3_lo = 1 - t3; w3_hi = t3;

        % Evaluate at corners and interpolate incrementally
        % Start with corner (lo, lo, lo)
        vals = contract_tucker(core, F1_lo, F2_lo, F3_lo, N, R) .* w1_lo .* w2_lo .* w3_lo;
        % Add corner (hi, lo, lo)
        vals = vals + contract_tucker(core, F1_hi, F2_lo, F3_lo, N, R) .* w1_hi .* w2_lo .* w3_lo;
        % Add corner (lo, hi, lo)
        vals = vals + contract_tucker(core, F1_lo, F2_hi, F3_lo, N, R) .* w1_lo .* w2_hi .* w3_lo;
        % Add corner (hi, hi, lo)
        vals = vals + contract_tucker(core, F1_hi, F2_hi, F3_lo, N, R) .* w1_hi .* w2_hi .* w3_lo;
        % Add corner (lo, lo, hi)
        vals = vals + contract_tucker(core, F1_lo, F2_lo, F3_hi, N, R) .* w1_lo .* w2_lo .* w3_hi;
        % Add corner (hi, lo, hi)
        vals = vals + contract_tucker(core, F1_hi, F2_lo, F3_hi, N, R) .* w1_hi .* w2_lo .* w3_hi;
        % Add corner (lo, hi, hi)
        vals = vals + contract_tucker(core, F1_lo, F2_hi, F3_hi, N, R) .* w1_lo .* w2_hi .* w3_hi;
        % Add corner (hi, hi, hi)
        vals = vals + contract_tucker(core, F1_hi, F2_hi, F3_hi, N, R) .* w1_hi .* w2_hi .* w3_hi;
    else
        error('Unknown method: %s. Use ''nearest'' or ''linear''.', method);
    end

    % Reshape output
    if length(sz) >= 3
        vals = reshape(vals, sz(1:end-1));
    end
end

function vals = contract_tucker(core, F1, F2, F3, N, R)
    % Contract Tucker core with factor matrices
    % vals(p) = sum_{i,j,k} core(i,j,k) * F1(p,i) * F2(p,j) * F3(p,k)

    % Step 1: Contract core with F3 over mode 3
    temp1 = tensorprod(core, F3, 3, 2);  % r1 x r2 x N

    % Step 2: Permute to N x r1 x r2
    temp1 = permute(temp1, [3, 1, 2]);  % N x r1 x r2

    % Contract over r1: sum_i temp1(p,i,j) * F1(p,i)
    temp2 = squeeze(sum(temp1 .* reshape(F1, N, R(1), 1), 2));  % N x r2

    % Contract over r2: sum_j temp2(p,j) * F2(p,j)
    vals = sum(temp2 .* F2, 2);  % N x 1
end

function F = qtt_batch_eval(qtt, subs)
    % Batch-evaluate a QTT tensor at multiple multi-indices
    % Input:
    %   qtt  - tt_tensor of dimension d+1, where the last mode has size r
    %          (joint QTT format: modes 1..d are binary, mode d+1 is the rank)
    %   subs - N x d matrix of multi-indices, each row is (i1,...,id) with i_k in {1,2}
    % Output:
    %   F    - N x r matrix, F(p,:) = qtt(subs(p,1), ..., subs(p,d), :)

    N = size(subs, 1);
    cores = core(qtt);
    r = qtt.r;  % TT ranks: r(1)=1, r(k) is rank between core k-1 and k, r(end)=1
    num_cores = length(cores);
    d = num_cores - 1;  % last core is the rank mode

    % TT-toolbox convention: core k has size n_k x r(k) x r(k+1)
    % First core: n1 x r(2) (2D, r(1)=1 implicit)
    % Intermediate: n_k x r(k) x r(k+1) (3D)
    % Last core: n_d x r(d) (2D, r(d+1)=1 implicit)

    % Pre-allocate state buffer with maximum rank
    max_rank = max(r);
    state = zeros(N, max_rank);

    % Initialize: start with first core sliced at subs(:,1)
    state(:, 1:r(2)) = cores{1}(subs(:,1), :);  % N x r(2)

    % Process intermediate binary cores (modes 2 to d)
    for k = 2:d
        % Select based on subs(:,k) and compute only needed contractions
        mask1 = (subs(:, k) == 1);

        % In-place update: compute contractions and store in state
        state(mask1, 1:r(k+1)) = state(mask1, 1:r(k)) * reshape(cores{k}(1, :, :), [r(k), r(k+1)]);
        state(~mask1, 1:r(k+1)) = state(~mask1, 1:r(k)) * reshape(cores{k}(2, :, :), [r(k), r(k+1)]);
    end

    % Process last core (rank mode)
    % F(p, j) = sum_i state(p, i) * Gd(j, i)
    F = state(:, 1:r(d+1)) * cores{num_cores}';  % N x Tucker_rank
end
