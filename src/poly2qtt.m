function [qtt,truevals] = poly2qtt(polycoef, varargin)
%POLY2QTT Convert polynomial(s) to QTT format.
%   QTT = POLY2QTT(POLYCOEF) constructs a QTT tensor representation of the
%   polynomial(s) specified by POLYCOEF on a uniform grid of 2^d points.
%
%   POLYCOEF can be:
%     - A matrix of size (m+1) x M, where each column contains coefficients
%       of a polynomial of degree m (in ascending degree order)
%     - A cell {matrix, type} where type is 'mono' or 'cheb' indicating
%       the polynomial basis (monomial or Chebyshev)
%
%   Optional arguments:
%     D         - QTT dimension (grid size 2^D), default: 20
%     DOM       - Domain [a,b], default: [-1,1]
%     'method'  - Approximation method:
%                 'ConstructiveLowRank' (default), 'ConstructiveFullRank',
%                 'ConstructiveMonomial', 'TTSVD', 'TTCross', 'Horner'
%     'tol'     - Truncation tolerance, default: 1e-12
%     'order'   - Position of polynomial index mode for multiple polynomials:
%                 'first' or 'last' (default)
%
%   [QTT, TRUEVALS] = POLY2QTT(...) also returns the true values on the grid
%   (only computed for TTSVD method).
%
%   See also: TT_TENSOR

%% --- Input Parsing ---
% Validation functions
validateD = @(d) isnumeric(d) && isscalar(d) && (d > 1) && (round(d) == d);
validateDom = @(dom) isnumeric(dom) && isvector(dom) && (length(dom) == 2) && (dom(1) <= dom(2));

% seperate positional arguments and name-value pairs
nv_pairs = {};
first_nv_pair_idx = find(cellfun(@(x) ischar(x) || isstring(x), varargin), 1);

if isempty(first_nv_pair_idx)
    positional_args = varargin;
else
    positional_args = varargin(1:first_nv_pair_idx-1);
    nv_pairs = varargin(first_nv_pair_idx:end);
end

defaultD = 20;
defaultDom = [-1, 1];

num_pos = length(positional_args);
if num_pos == 2
    d = positional_args{1};
    dom = positional_args{2};
elseif num_pos == 0
    d = defaultD;
    dom = defaultDom;
elseif num_pos > 2
    error('poly2qtt:TooManyPositionalArgs', 'Too many positional arguments. Expected at most D and DOMAIN.');
elseif num_pos == 1
    arg1 = positional_args{1};
    if validateD(arg1)
        d = arg1;
        dom = defaultDom;
    elseif validateDom(arg1)
        dom = arg1;
        d = defaultD;
    else
        error('poly2qtt:InvalidFirstOptionalArg', 'The first optional variable must be D or DOMAIN.');
    end
end


% Default values
defaultQ = 2;
defaultTol = 1e-12;
defaultDebug = false;
defaultMethod = 'ConstructiveLowRank';
expectedMethods = {'TTSVD', 'TTCross', 'Horner', ...
    'ConstructiveLowRank', 'ConstructiveFullRank', 'ConstructiveMonomial'};
defaultOrder = 'last';
expectedOrders = {'first', 'last'};
% If order is 'first', then the dimension of the polynomial is added as the first dimension of the QTT
% i.e., the QTT will be of size M x 2 x 2 x ... x 2 (d times)
% If order is 'last', then the QTT will be of size 2 x 2 x ... x 2 (d times) x M

% Validation functions for nv pairs
% polycoef: either a matrix of coefficients or {matrix, type} with type == 'mono' or 'cheb'
% matrix of coefficients: size (m+1) x M, where each column is the coefficients (in ascending degree order) of a polynomial of degree m
validatePoly = @(f) ismatrix(f) || (iscell(f) && ~isempty(f) && ismatrix(f{1}) && (isscalar(f) || (length(f) == 2 && any(strcmpi(f{2}, {'mono', 'cheb'})))));
% quantics: either a scalar or a vector of length d, all integers > 1
validateQ = @(q) isnumeric(q) && (isscalar(q) || (isvector(q) && length(q) == d)) && all(q > 1) && all(round(q) == q);
% tol: positive scalar
validateTol = @(x) isnumeric(x) && isscalar(x) && (x > 0);
% debug: logical or 0 or 1
validateDebug = @(x) islogical(x) || (isnumeric(x) && isscalar(x) && (x == 0 || x == 1));

% Setup inputParser for 'quantics'
p = inputParser;
p.FunctionName = 'poly2qtt';

% Required 'polycoef' argument
addRequired(p, 'polycoef', validatePoly);
% Optional positional arguments
addOptional(p, 'd', defaultD, validateD);
addOptional(p, 'dom', defaultDom, validateDom);
% Name-Value pairs
addParameter(p, 'quantics', defaultQ, validateQ);
addParameter(p, 'tol', defaultTol, validateTol);
addParameter(p, 'debug', defaultDebug, validateDebug);
addParameter(p, 'method', defaultMethod, @(x) any(validatestring(x,expectedMethods)));
addParameter(p, 'order', defaultOrder, @(x) any(validatestring(x,expectedOrders)));

% Parse the inputs (polycoef + remaining varargin for 'quantics')
parse(p, polycoef, d, dom, nv_pairs{:});

% Retrieve the final parsed values
polycoef_input = p.Results.polycoef;
d = p.Results.d;
dom = p.Results.dom;
q = p.Results.quantics;
tol = p.Results.tol;
debug = p.Results.debug;
method = p.Results.method;
order = p.Results.order;

polytype = 'cheb'; % default type
% Parse polycoef to extract coefficients matrix and polynomial type
if iscell(polycoef_input)
    polycoef = polycoef_input{1};
    if length(polycoef_input) >= 2
        polytype = polycoef_input{2};
    end
else
    polycoef = polycoef_input;
end
polyfun = poly2handle(polycoef, polytype, dom);

M = size(polycoef, 2); % number of polynomials
l = size(polycoef, 1); % polynomial degree + 1
m = l - 1; % polynomial degree

% Not implemented yet for mode size other than 2
if length(q) > 1 || q > 2
    q = 2;
    warning('poly2qtt:NotImplemented', 'Currently, only quantics of size 2 are implemented. Using default quantics = 2.');
end

last = ifelse(strcmpi(order, 'first'), false, true);

qtt_size = q*ones(1,d);
qtt_dim = d;
if M > 1
    qtt_dim = d + 1;
    qtt_size = ifelse(last, [qtt_size, M], [M, qtt_size]);
end
% --- End of Input Parsing ---

%% --- Main function logic starts here ---

%% Corner cases: m == 0
if m == 0
    r = ones(qtt_dim+1,1);
    if M == 1
        cr = [polycoef*ones(2,1);ones(qtt_dim*2-2,1)];
    else
        Mn = polycoef .* ones(2,1); % size M x 2
        [Mn, R] = qr(Mn', 'econ');
        if last
            r(qtt_dim) = 2;
            Mn = Mn'; R = R';
            cr = [ones((qtt_dim-2)*2,1);R(:);Mn(:)];
        else
            r(2) = 2;
            cr = [Mn(:);R(:);ones((qtt_dim-2)*2,1)];
        end
    end
    truevals = NaN;
    qtt = tt_tensor();
    qtt.d    = qtt_dim;
    qtt.r    = r;
    qtt.n    = qtt_size;
    qtt.core = cr;
    qtt.ps   = cumsum([1;qtt_size'.*r(1:qtt_dim).*r(2:qtt_dim+1)]);
    return
end

N = 2^d;
h = (dom(2) - dom(1)) / N;
truevals = NaN; % Not computed by default

%% Use TTSVD or TTCross methods
if strcmpi(method, 'TTSVD')
    % Use TTSVD method — polyeval avoids O(N*(m+1)) intermediate
    grid = linspace(dom(1) + h/2, dom(2) - h/2, N)';
    truevals = polyeval(grid, polycoef, polytype, dom); % size N x M
    qtt = tt_tensor(reshape(ifelse(last,truevals,truevals'),qtt_size), tol);
    return;
elseif strcmpi(method, 'TTCross')
    % Use TTCross method
    qtt = dmrg_cross(qtt_dim,qtt_size,@(ind) fun2tensor(polyfun,ind,dom,d,order,M),tol,'vec',true,'verb',false);
    return;
end

%% Use constructive methods
% --- Define the quantics grid point x^(k)_ik ---
% First define the constant offset
% sum(X_offset) should equals -1 + h/2
X_offset = zeros(d,1);
X_offset(d) = (dom(1) + h/2)/2;
X_offset(2:(d-1)) = (dom(1) + h/2)/2/(d-2);
% size(X) = dxq
% X(k,i_k) = x^(k)_ik for k = 1,...,d and i_k = 1,...,q
X = [X_offset, X_offset + q.^(0:(d-1))' * h];
if last
    ind_list = [qtt_dim,d,(d-1):-1:2,1];
    X = flipud(X);
else % order == 'first'
    ind_list = [1,2,(qtt_dim+2-d):qtt_dim];
end
Igeq = cumsum(X,'reverse');

% --- Initialize the cores cell array ---
cores = cell(qtt_dim,1);

if strcmpi(method, 'ConstructiveMonomial')
    % This is only for experimental purposes and should never be used
    % in practice since it is extremely unstable.
    % It only takes the input as a polynomial in monomial basis
    if ~strcmpi(polytype, 'mono')
        error('poly2qtt:InvalidInput', 'Input must be a polynomial in monomial basis.');
    end
    % --- extra cores due to multiple polynomials ---
    % if last, they are cores{d+1} and cores{d}
    % if first, they are cores{1} and cores{2}
    Mm = zeros(l, l, M);
    for j = 1:M
        Mm(:, :, j) = hankel(polycoef(:,j));
    end
    Mm = Mm .* pascal(l);
    if last; Mm = permute(Mm, [2,1,3]); end
    Mm = tensorprod(X(1,:)'.^((0:m)), Mm, 2, 1); % size 2 x l x M
    if M == 1
        % If single polynomial
        % update cores{1} (first) or cores{d} (last)
        cores{ind_list(1)} = Mm; % size 2 x l
    else
        Mm = reshape(Mm, [2*l, M])';
        [Mm, R] = qr(Mm, 'econ');
        cores{ind_list(1)} = Mm; % size M x r with r = min(M,2*l)
        cores{ind_list(2)} = permute(reshape(R, [min(M,2*l), 2, l]), [2, 1, 3]); % size 2 x r x l
    end

    % --- Other cores ---
    % if last, they are cores{1} to cores{d-1}
    % if first, they are cores{3} to cores{d+1}
    P = abs(pascal(l,1));
    corek = zeros(2,l,l);
    for k = 3:d
        corek(1,:,:) = flipud(hankel(X(k-1,1).^(m:-1:0))).* P;
        corek(2,:,:) = flipud(hankel(X(k-1,2).^(m:-1:0))).* P;
        cores{ind_list(k)} = corek;
    end
    cores{ind_list(d+1)} = X(d,:)'.^(0:m);

elseif strcmpi(method, 'ConstructiveFullRank')
    % ConstructiveFullRank
    CC = chebvals2coeffsMat(l);
    % --- extra cores due to multiple polynomials ---
    % if last, they are cores{d+1} and cores{d}
    % if first, they are cores{1} and cores{2}
    xx = X(1,:)'; yy = chebpts(l,Igeq(2,:));
    Mc = polyfun(xx+yy'); % 2 x l x M
    Mc = tensorprod(Mc,CC,2,2); % = but faster than M * CC.';
    if M == 1
        % If single polynomial
        % update cores{1} (first) or cores{d} (last)
        cores{ind_list(1)} = Mc; % size 2 x l
    else
        Mc = reshape(permute(Mc,[1,3,2]), [2*l, M])';
        [Mc, R] = qr(Mc, 'econ');
        cores{ind_list(1)} = Mc; % size M x r with r = min(M,2*l)
        cores{ind_list(2)} = permute(reshape(R, [min(M,2*l), 2, l]), [2, 1, 3]); % size 2 x r x l
    end

    % --- Other cores ---
    % if last, they are cores{1} to cores{d-1}
    % if first, they are cores{3} to cores{d+1}
    for k = 3:d
        xx = X(k-1,:)'; yy = chebpts(l,Igeq(k,:));
        C = chebpolyeval(xx+yy',reshape(0:m,[1,1,l]),Igeq(k-1,:));
        cores{ind_list(k)} = tensorprod(C,CC,2,2);
    end
    cores{ind_list(d+1)} = squeeze(chebpolyeval(X(d,:)',reshape(0:m,[1,1,l]),Igeq(d,:)));

elseif strcmpi(method, 'ConstructiveLowRank')
    % ConstructiveLowRank
    tol = tol/((d-1)*(sqrt(l) + 1));
    r = [1;zeros(qtt_dim-1,1);1];
    CC = chebvals2coeffsMat(l);
    xx = X(1,:)'; yy = chebpts(l,Igeq(2,:));
    Mc = polyfun(xx+yy'); % 2 x l x M
    if M == 1
        % If single polynomial
        % update cores{1} (first) or cores{d} (last)
        [U,S,V, r(ind_list(3))] = truncated_svd(Mc, tol);
        V = CC * V * S;
        cores{ind_list(1)} = U; % size 2 x r1
    else
        Mc = reshape(Mc,[2*l, M])'; % M x l*2
        [U,S,V,r(ind_list(2))] = truncated_svd(Mc, tol);
        V = V * S; % 2*l x r2
        cores{ind_list(1)} = U; % size M x r2

        V = reshape(permute(reshape(V, [2, l, r(ind_list(2))]), [3,1,2]), [r(ind_list(2))*2, l]); % size r2*2 x l
        [U,S,V,r(ind_list(3))] = truncated_svd(V, tol);
        V = CC * V * S;
        cores{ind_list(2)} = permute(reshape(U, [r(ind_list(2)),2,r(ind_list(3))]),[2,1,3]); % size 2 x r x r1
    end

    for k = 3:d
        xx = X(k-1,:)'; yy = chebpts(l,Igeq(k,:));
        xall = xx + yy'; % 2 x l
        C = reshape(clenshaw_eval_mex(xall(:), V, Igeq(k-1,1), Igeq(k-1,2)), [2, l, size(V,2)]); % 2 x l x r1
        C = reshape(permute(C,[3,1,2]), [2*r(ind_list(k)), l]);

        [U,S,V, r(ind_list(k+1))] = truncated_svd(C, tol);
        V = CC * V * S;
        cores{ind_list(k)} = permute(reshape(U, [r(ind_list(k)),2,r(ind_list(k+1))]),[2,1,3]); % size 2 x r1 x r2
    end
    cores{ind_list(d+1)} = clenshaw_eval_mex(X(d,:)', V, Igeq(d,1), Igeq(d,2)); % size 2 x r1
elseif strcmpi(method, 'Horner')
    % Horner's method
    x_qtt = poly2qtt({[zeros(1,M);ones(1,M)],'mono'}, d, dom, 'method', 'ConstructiveMonomial', 'order', order);
    if strcmpi(polytype, 'mono')
        bm = poly2qtt({polycoef(m:m+1,:),'mono'}, d, dom, 'method', 'ConstructiveMonomial', 'order', order);
        for ii = m-1:-1:1
            bm_1 = poly2qtt(polycoef(ii,:), d, dom, 'order', order) + x_qtt .* bm;
            bm = round(bm_1, tol);
        end
    elseif strcmpi(polytype, 'cheb')
        bmp1 = 0;
        bm = poly2qtt(polycoef(m+1,:), d, dom, 'order', order);
        for ii = m:-1:2
            bm_1 = poly2qtt(polycoef(ii,:), d, dom, 'order', order) + 2 * x_qtt .* bm - bmp1;
            bmp1 = bm;
            bm = round(bm_1, tol);
        end
        bm = poly2qtt(polycoef(1,:), d, dom, 'order', order) + x_qtt .* bm - bmp1;
    end
    qtt = round(bm, tol);
    return
end

if last
    for k = 2:(qtt_dim-1); cores{k} = permute(cores{k}, [1,3,2]); end
end

if strcmpi(method, 'ConstructiveLowRank')
    % This final rounding step calls custom tt_round function
    % to avoid the QR orthogonalization step
    if last; r = [1;r(1:qtt_dim-1);1]; end
    % The tol here is basically the user tol (default = 1e-12) divided by sqrt(d-1).
    % Theoretically we should use the same tol as before
    % but in practice we find that this leads to unnecessarily high ranks.
    % Use a higher tol here seems to be sufficient to maintain the overall error.
    cores = tt_round(cores, tol*sqrt(d-1)*(sqrt(l) + 1), ifelse(last,'lr','rl'), qtt_size, r);
end
%% --- Finalize the QTT tensor ---
qtt = tt_tensor(cores);
end

function polyfun = poly2handle(polycoef, polytype, dom)
    if strcmpi(polytype, 'mono')
        polyfun = @(x) reshape(horner_eval_mex(x(:), polycoef), [size(x), size(polycoef,2)]);
    elseif strcmpi(polytype, 'cheb')
        polyfun = @(x) reshape(clenshaw_eval_mex(x(:), polycoef, dom(1), dom(2)), [size(x), size(polycoef,2)]);
    end
end

function val = ifelse(cond, valTrue, valFalse)
    if cond
        val = valTrue;
    else
        val = valFalse;
    end
end

function [U, S, V, r] = truncated_svd(X, tol)
    [U, S, V] = svd(X, 'econ');
    r = sum(diag(S)./S(1,1) >= tol);
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);
end

function cores = tt_round(cores, tol, direction, n, r)
    % TT rounding to reduce ranks of TT cores
    % This does not have the QR orthogonalization step as in
    % the standard TT rounding algorithm since the cores are already
    % left- or right- orthogonalized in the constructive algorithms.
    % direction: 'lr' for left-to-right, 'rl' for right-to-left
    d = numel(cores);

    if nargin < 5
        n = zeros(1, d);
        r = ones(1, d + 1);
        for i = 1:d
            s = size(cores{i});
            n(i) = s(1);
            if numel(s) > 2
                r(i+1) = s(3);
            end
            if i > 1
                r(i) = s(2);
            else
                if numel(s) > 1
                    r(i+1) = s(2);
                end
            end
        end
    end

    if strcmpi(direction, 'lr')
        for ii = 1:(d-1)
            core0 = reshape(cores{ii}, [n(ii)*r(ii), r(ii+1)]);
            [U,S,V,r(ii+1)] = truncated_svd(core0, tol);
            cores{ii} = squeeze(reshape(U, [n(ii), r(ii), r(ii+1)])); % size n x r(i) x r(i+1)
            cores{ii+1} = permute(tensorprod(S*V', cores{ii+1}, 2, 2),[2,1,3]);
        end
    elseif strcmpi(direction, 'rl')
        for ii = d:-1:2
            core0 = reshape(permute(cores{ii}, [2,3,1]), [r(ii), r(ii+1)*n(ii)]);
            [U,S,V,r(ii)] = truncated_svd(core0, tol);
            cores{ii} = permute(reshape(V', [r(ii), r(ii+1), n(ii)]), [3,1,2]); % size n x r(i) x r(i+1)
            cores{ii-1} = tensorprod(cores{ii-1}, U*S, ndims(cores{ii-1}), 1);
        end
    else
        error('tt_round:InvalidDirection', 'Direction must be ''lr'' or ''rl''.');
    end
end


function y = chebpolyeval(x, k, I)
%CHEBPOLYEVAL Evaluate Chebyshev polynomials at given points.
%   Y = CHEBPOLYEVAL(X, K, I) evaluates the Chebyshev polynomial T_k(x)
%   at the points X, where the polynomials are defined on the interval I.

y = cos(k.*acos(I2unitI(x,I)));

end