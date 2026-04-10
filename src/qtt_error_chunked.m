function err = qtt_error_chunked(tt, polycoef, polytype, dom, varargin)
%QTT_ERROR_CHUNKED Compute inf-norm error between QTT(s) and polynomial.
%   ERR = QTT_ERROR_CHUNKED(TT, POLYCOEF, POLYTYPE, DOM) computes the
%   infinity-norm error ||full(TT) - polyeval(grid, POLYCOEF, POLYTYPE, DOM)||_inf
%   using chunked evaluation to avoid materializing the full 2^d grid.
%
%   TT can be a single tt_tensor or a cell array of tt_tensors.
%   When TT is a cell array, the ground truth is evaluated only once
%   per chunk and ERR is returned as a vector of length numel(TT).
%
%   ERR = QTT_ERROR_CHUNKED(..., 'chunkbits', C) sets the chunk size
%   to 2^C grid entries. Default: min(d, 25) (~256 MB per chunk).
%
%   ERR = QTT_ERROR_CHUNKED(..., 'order', 'last'|'first') specifies
%   the polynomial index mode position for M>1. Default: 'last'.
%
%   ERR = QTT_ERROR_CHUNKED(..., 'verbose', true) prints progress.
%
%   This function is useful for large d where full(TT) would exhaust
%   available memory. For small d (e.g., d <= 25), direct error
%   computation via full(TT) is usually faster.
%
%   Example:
%     tt = poly2qtt({coefs, 'cheb'}, 30, 'method', 'ConstructiveLowRank');
%     err = qtt_error_chunked(tt, coefs, 'cheb', [-1,1]);
%
%   See also: QTT_CHUNK_EVAL, POLYEVAL, POLY2QTT

% Parse inputs
p = inputParser;
addParameter(p, 'chunkbits', [], @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'order', 'last', @(x) any(validatestring(x, {'first', 'last'})));
addParameter(p, 'verbose', false, @islogical);
parse(p, varargin{:});
chunk_bits = p.Results.chunkbits;
order = validatestring(p.Results.order, {'first', 'last'});
verbose = p.Results.verbose;

% Handle single TT or cell array of TTs
if ~iscell(tt)
    tt = {tt};
    was_scalar = true;
else
    was_scalar = false;
end
n_tt = numel(tt);

% Determine M and d
M = size(polycoef, 2);
d_tt = tt{1}.d;
if M == 1
    d = d_tt;
else
    d = d_tt - 1;
end

if isempty(chunk_bits)
    chunk_bits = min(d, 25);
end
c = chunk_bits;

N = 2^d;
h = (dom(2) - dom(1)) / N;
n_chunks = 2^(d - c);

err = zeros(n_tt, 1);
t_start = tic;

for k = 0:(n_chunks - 1)
    % Compute ground truth for this chunk (once for all TTs)
    idx_start = k * 2^c;
    grid_chunk = dom(1) + h/2 + h * (idx_start:(idx_start + 2^c - 1))';
    truth_chunk = polyeval(grid_chunk, polycoef, polytype, dom);

    % Compare each TT against the same truth
    for j = 1:n_tt
        qtt_vals = qtt_chunk_eval(tt{j}, k, 'chunkbits', c, 'order', order, 'M', M);
        err(j) = max(err(j), max(abs(qtt_vals(:) - truth_chunk(:))));
    end

    if verbose && (mod(k+1, max(1, floor(n_chunks/10))) == 0 || k == n_chunks - 1)
        fprintf('  Chunk %d/%d (%.1f%%), elapsed %.1fs, current max err: [%s]\n', ...
            k+1, n_chunks, (k+1)/n_chunks*100, toc(t_start), ...
            strjoin(arrayfun(@(e) sprintf('%.4e',e), err, 'UniformOutput', false), ', '));
    end
end

if was_scalar
    err = err(1);
end

end
