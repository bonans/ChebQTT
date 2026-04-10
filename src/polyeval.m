function y = polyeval(x, V, polytype, dom)
%POLYEVAL Evaluate polynomial linear combination at given points.
%   Y = POLYEVAL(X, V) evaluates monomial combination at points X using
%   Horner's scheme (MEX). V is (m+1) x M coefficient matrix.
%
%   Y = POLYEVAL(X, V, POLYTYPE) specifies the basis: 'mono' (default)
%   or 'cheb'.
%
%   Y = POLYEVAL(X, V, POLYTYPE, DOM) specifies the domain [a,b] for
%   Chebyshev polynomials (default: [-1,1]).
%
%   Uses compiled MEX routines (Horner / Clenshaw) for all evaluations.
%   These are consistently faster than vectorized or native MATLAB across
%   all tested configurations (N, m, M).
%
%   See also: POLY2QTT, HORNER_EVAL_MEX, CLENSHAW_EVAL_MEX

if nargin < 3 || isempty(polytype)
    polytype = 'mono';
end
if nargin < 4 || isempty(dom)
    dom = [-1, 1];
end

% Determine input shape
orig_size = size(x);
d = ndims(x);
M = size(V, 2);
is_vec = (d == 2 && isvector(x));
if is_vec
    d = 1;
    x = x(:);
end

if strcmpi(polytype, 'mono')
    y = horner_eval_mex(x(:), V);
elseif strcmpi(polytype, 'cheb')
    y = clenshaw_eval_mex(x(:), V, dom(1), dom(2));
else
    error('polyeval:InvalidType', 'polytype must be ''mono'' or ''cheb''.');
end

% Reshape output for N-d input arrays
if M > 1 && ~is_vec && d > 1
    y = reshape(y, [orig_size, M]);
elseif M == 1 && is_vec
    y = reshape(y, orig_size);
end

end
