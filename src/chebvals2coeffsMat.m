function F = chebvals2coeffsMat(n)
%CHEBVALS2COEFFSMAT Chebyshev values to coefficients transformation matrix.
%   F = CHEBVALS2COEFFSMAT(N) returns the transformation matrix that maps
%   function evaluations at N Chebyshev nodes to Chebyshev coefficients.
%
%   See also: CHEBCOMBEVAL

if n < 2 || floor(n) ~= n
    error('chebvals2coeffsMat:InvalidInput', 'n must be an integer >= 2');
end
F = (2/(n-1)).*cos( (0:n-1)' * ((n-1:-1:0).*pi/(n-1)) );

F(:,[1,n]) = F(:,[1,n])/2;
F([1,n],:) = F([1,n],:)/2;
end
