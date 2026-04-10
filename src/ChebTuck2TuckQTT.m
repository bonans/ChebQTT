function TuckQTT = ChebTuck2TuckQTT(ChebTuck, varargin)
%CHEBTUCK2TUCKQTT Convert a Chebfun3 object to QTT-Tucker format.
%   TUCKQTT = CHEBTUCK2TUCKQTT(CHEBTUCK) converts a Chebfun3 object in
%   ChebTuck format to a fully discrete QTT-Tucker-like format.
%
%   Optional name-value pairs:
%     'tol'    - Tolerance for QTT approximation (default: 1e-12)
%     'd'      - Number of quantics dimensions (default: 20)
%     'order'  - 'first' or 'last' for polynomial index mode (default: 'last')
%     'method' - Method for QTT approximation (default: 'ConstructiveLowRank')
%
%   Output:
%     TUCKQTT - Struct with fields:
%       .core  - Core tensor in Tucker format (full format)
%       .cols  - Joint QTT format for polynomials in first mode
%       .rows  - Joint QTT format for polynomials in second mode
%       .tubes - Joint QTT format for polynomials in third mode
%       .dom   - Domain of the function (3x2 matrix)
%       .d     - Number of quantics dimensions (3x1 vector)
%
%   See also: POLY2QTT, TUCKQTT_EVAL

    % Input validation
    if ~isa(ChebTuck, 'chebfun3')
        error('ChebTuck2TuckQTT:InvalidInput', 'Input must be a chebfun3 object.');
    end

    if nargin < 2
        varargin = {};
    else
        varargin = varargin(:)';
    end

    % Default value for d
    d = 20;
    for k = 1:2:length(varargin)
        key = varargin{k};
        value = varargin{k+1};
        if strcmpi(key, 'd')
            d = value;
        end
    end
    TuckQTT.d = d * ones(3,1);
    cols_polycoef = {chebcoeffs(ChebTuck.cols), 'cheb'};
    rows_polycoef = {chebcoeffs(ChebTuck.rows), 'cheb'};
    tubes_polycoef = {chebcoeffs(ChebTuck.tubes), 'cheb'};
    TuckQTT.dom = reshape(ChebTuck.domain, [2,3])';
    TuckQTT.core = ChebTuck.core;
    TuckQTT.cols = poly2qtt(cols_polycoef, varargin{:});
    TuckQTT.rows = poly2qtt(rows_polycoef, varargin{:});
    TuckQTT.tubes = poly2qtt(tubes_polycoef, varargin{:});
end
