function vals = qtt_chunk_eval(tt, chunk_idx, varargin)
%QTT_CHUNK_EVAL Extract a chunk of entries from a QTT tensor.
%   VALS = QTT_CHUNK_EVAL(TT, CHUNK_IDX) extracts 2^c contiguous grid
%   entries from the QTT tensor TT at chunk index CHUNK_IDX (0-based).
%   The chunk size 2^c defaults to min(d, 25) grid bits.
%
%   VALS = QTT_CHUNK_EVAL(TT, CHUNK_IDX, 'chunkbits', C) specifies the
%   number of grid bits per chunk. Each chunk contains 2^C entries.
%
%   VALS = QTT_CHUNK_EVAL(TT, CHUNK_IDX, 'order', 'last') specifies the
%   position of the polynomial index mode for multi-polynomial QTTs.
%   Default: 'last'. Only relevant when the QTT has M > 1.
%
%   VALS = QTT_CHUNK_EVAL(TT, CHUNK_IDX, 'M', M) specifies the number of
%   polynomials, allowing automatic detection of grid dimension d = tt.d
%   (M=1) or d = tt.d - 1 (M>1). Default: 1.
%
%   The output VALS has size [2^C, M] where M is the number of polynomials.
%
%   Grid ordering: chunk CHUNK_IDX covers grid indices
%   CHUNK_IDX*2^C + 1 through (CHUNK_IDX+1)*2^C (1-based), matching the
%   ordering of full(TT).
%
%   Example (M=1):
%     % Full evaluation via chunking for d=30
%     d = 30; c = 25; n_chunks = 2^(d - c);
%     vals = zeros(2^d, 1);
%     for k = 0:(n_chunks-1)
%         vals(k*2^c+1 : (k+1)*2^c) = qtt_chunk_eval(tt, k, 'chunkbits', c);
%     end
%     % vals should match full(tt) (up to reshaping)
%
%   See also: POLY2QTT, POLYEVAL, TT_TENSOR

% Parse inputs
p = inputParser;
addParameter(p, 'chunkbits', [], @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'order', 'last', @(x) any(validatestring(x, {'first', 'last'})));
addParameter(p, 'M', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parse(p, varargin{:});
order = validatestring(p.Results.order, {'first', 'last'});
M = p.Results.M;
chunk_bits = p.Results.chunkbits;

% Get TT info
d_tt = tt.d;
cores = core2cell(tt); % each cell: r(k) × n(k) × r(k+1)

% Determine grid dimension
if M == 1
    d = d_tt;
else
    d = d_tt - 1;
end

if isempty(chunk_bits)
    chunk_bits = min(d, 25); % default: ~256 MB per chunk
end
c = chunk_bits;

assert(c >= 1 && c <= d, 'qtt_chunk_eval: chunkbits must satisfy 1 <= chunkbits <= d (=%d)', d);
assert(chunk_idx >= 0 && chunk_idx < 2^(d-c), ...
    'qtt_chunk_eval: chunk_idx must satisfy 0 <= chunk_idx < 2^(d-c) (=%d)', 2^(d-c));

%% Dispatch based on M and order
if M == 1
    vals = extract_M1(cores, d, c, chunk_idx);
elseif strcmpi(order, 'last')
    vals = extract_joint_last(cores, d, c, chunk_idx, M);
else % order == 'first'
    vals = extract_joint_first(cores, d, c, chunk_idx, M);
end

end

%% ---- M=1: all cores are grid modes ----
function vals = extract_M1(cores, d, c, chunk_idx)
    if c == d
        % No upper bits: expand all cores
        result = expand_left(cores, 1, d);
        vals = result(:); % [2^d, 1]
        return;
    end

    % Right contraction: cores d → c+1 at fixed bit indices
    v = contract_right(cores, c+1, d, chunk_idx, c);

    % Left expansion: cores 1 → c
    result = expand_left(cores, 1, c);

    % Final contraction
    vals = result * v;
end

%% ---- M>1, order='last': grid = cores 1:d, poly = core d+1 ----
function vals = extract_joint_last(cores, d, c, chunk_idx, M)
    d_tt = d + 1;

    % Start from poly core (rightmost): r(d+1) × M × 1
    V = squeeze(cores{d_tt}); % r(d+1) × M
    if M == 1; V = V(:); end  % ensure column vector

    if c == d
        % No upper grid bits: expand all grid cores and multiply by poly core
        result = expand_left(cores, 1, d);
        vals = result * V; % [2^d, M]
        return;
    end

    % Contract upper grid cores d → c+1 at fixed bits
    bits = bitget(chunk_idx, 1:(d-c)) + 1; % 1-indexed
    for j = d:-1:(c+1)
        G = squeeze(cores{j}(:, bits(j-c), :)); % r(j) × r(j+1)
        if isscalar(G); G = reshape(G, [1, 1]); end
        V = G * V; % r(j) × M
    end
    % V is now r(c+1) × M

    % Left expansion: cores 1 → c
    result = expand_left(cores, 1, c);
    % result is 2^c × r(c+1)

    % Final: 2^c × M
    vals = result * V;
end

%% ---- M>1, order='first': poly = core 1, grid = cores 2:d+1 ----
function vals = extract_joint_first(cores, d, c, chunk_idx, M)
    d_tt = d + 1;

    if c == d
        % No upper grid bits: expand all cores (poly + grid)
        result = expand_left(cores, 1, d_tt);
        vals = reshape(result, [M, 2^d])';
        return;
    end

    % Right contraction: upper grid cores (d+1) → (c+2) at fixed bits
    % Grid bits c+1 through d correspond to TT cores c+2 through d+1
    bits = bitget(chunk_idx, 1:(d-c)) + 1; % 1-indexed
    % Start from rightmost upper grid core = TT core d+1
    v = squeeze(cores{d_tt}(:, bits(d-c), :)); % r(d+1) × 1 (last TT core)
    v = v(:);
    for j = d:-1:(c+2)
        % TT core j: upper grid bit j-1, bit index = (j-1) - c = j-c-1
        G = squeeze(cores{j}(:, bits(j-c-1), :)); % r(j) × r(j+1)
        if isscalar(G); G = reshape(G, [1, 1]); end
        v = G * v; % r(j) × 1
    end
    % v is now r(c+2) × 1

    % Left expansion: poly core 1 + lower grid cores 2 → c+1
    result = expand_left(cores, 1, c+1);
    % result is M * 2^c × r(c+2)

    % Contract with right boundary
    result = result * v; % M * 2^c × 1

    % Reshape: column-major [M, 2^c] → transpose to [2^c, M]
    vals = reshape(result, [M, 2^c])';
end

%% ---- Helper: contract cores from right at fixed indices ----
function v = contract_right(cores, start_core, end_core, chunk_idx, c)
    % Contract TT cores end_core down to start_core at fixed bit indices
    % determined by chunk_idx.
    % start_core corresponds to grid bit c+1, end_core to grid bit d.
    d = end_core;   % last grid core index
    n_upper = d - c; % number of upper bits
    bits = bitget(chunk_idx, 1:n_upper) + 1; % 1-indexed

    % Start from rightmost core
    v = squeeze(cores{end_core}(:, bits(n_upper), :));
    v = v(:); % column vector

    for j = (end_core-1):-1:start_core
        bit_idx = j - c; % which bit in the bits array
        G = squeeze(cores{j}(:, bits(bit_idx), :));
        if isscalar(G); G = reshape(G, [1, 1]); end
        v = G * v;
    end
end

%% ---- Helper: expand cores from left ----
function result = expand_left(cores, start_core, end_core)
    % Expand TT cores start_core through end_core from left.
    % Mimics the full() algorithm for a sub-TT.
    % Returns matrix of size [prod(n(start_core:end_core)), r(end_core+1)]

    % Core start_core: r(start) × n × r(start+1)
    % For the first TT core (start=1), r(1) = 1
    r1 = size(cores{start_core}, 1);
    n1 = size(cores{start_core}, 2);
    r2 = size(cores{start_core}, 3);
    result = reshape(cores{start_core}, [r1 * n1, r2]);

    for j = (start_core+1):end_core
        rj = size(cores{j}, 1);
        nj = size(cores{j}, 2);
        rj1 = size(cores{j}, 3);
        cr = reshape(cores{j}, [rj, nj * rj1]);
        result = reshape(result, [numel(result) / rj, rj]);
        result = result * cr;
    end
    % After loop: result is [prod(n(start:end-1)), n(end)*r(end+1)]
    % Final reshape to separate the last n and r dimensions:
    % [r(start)*prod(n(start:end)), r(end+1)]
    r_end = size(cores{end_core}, 3); % r(end_core + 1)
    result = reshape(result, [], r_end);
end
