function err = qtt_poly_error(tt, polycoef, polytype, dom, varargin)
%QTT_POLY_ERROR Compute inf-norm error between QTT(s) and polynomial.
%   ERR = QTT_POLY_ERROR(TT, POLYCOEF, POLYTYPE, DOM) computes the
%   infinity-norm error between the QTT tensor TT and the polynomial
%   defined by POLYCOEF on the midpoint grid over DOM.
%
%   TT can be a single tt_tensor or a cell array of tt_tensors.
%   When TT is a cell array, all tensors must approximate the same
%   polynomial (same d). The ground truth is evaluated only once and
%   ERR is returned as a vector of length numel(TT).
%
%   Automatically selects direct full() evaluation for small d,
%   or chunked evaluation for large d.
%
%   ERR = QTT_POLY_ERROR(..., 'order', 'first'|'last') specifies
%   the polynomial index mode position for M>1. Default: 'last'.
%
%   ERR = QTT_POLY_ERROR(..., 'threshold', T) sets the d threshold
%   above which chunked evaluation is used. Default: 25.
%
%   ERR = QTT_POLY_ERROR(..., 'verbose', true) prints progress
%   (only relevant for chunked evaluation).
%
%   See also: QTT_ERROR_CHUNKED, POLYEVAL, POLY2QTT

p = inputParser;
addParameter(p, 'order', 'last', @(x) any(validatestring(x, {'first','last'})));
addParameter(p, 'threshold', 25, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'verbose', false, @islogical);
parse(p, varargin{:});
order     = validatestring(p.Results.order, {'first','last'});
threshold = p.Results.threshold;
verbose   = p.Results.verbose;

% Handle single TT or cell array of TTs
if ~iscell(tt)
    tt = {tt};
    was_scalar = true;
else
    was_scalar = false;
end
n_tt = numel(tt);

M = size(polycoef, 2);
d_tt = tt{1}.d;
if M == 1
    d = d_tt;
else
    d = d_tt - 1;
end

err = zeros(n_tt, 1);

if d <= threshold
    %% Direct evaluation via full() — ground truth computed once
    N = 2^d;
    h = (dom(2) - dom(1)) / N;
    grid_pts = dom(1) + h/2 + h*(0:N-1)';
    truevals = polyeval(grid_pts, polycoef, polytype, dom);

    for k = 1:n_tt
        tt_full = full(tt{k});
        if M == 1
            tt_vals = tt_full(:);
        elseif strcmpi(order, 'first')
            tt_vals = reshape(tt_full, [M, N])';
        else  % last
            tt_vals = reshape(tt_full, [N, M]);
        end
        err(k) = norm(tt_vals(:) - truevals(:), 'inf');
    end
else
    %% Chunked evaluation (memory-efficient for large d)
    err = qtt_error_chunked(tt, polycoef, polytype, dom, ...
        'order', order, 'verbose', verbose,'chunkbits', threshold);
end

if was_scalar
    err = err(1);
end

end
