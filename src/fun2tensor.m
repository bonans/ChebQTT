function fun_tensor = fun2tensor(fun, ind, dom, d, order, M)
%FUN2TENSOR Evaluate function(s) on a tensorized grid for TT-cross.
%   FUN_TENSOR = FUN2TENSOR(FUN, IND, DOM, D) evaluates a single function
%   at specific indices of the QTT tensor on a uniform grid.
%
%   FUN_TENSOR = FUN2TENSOR(FUN, IND, DOM, D, ORDER, M) handles multiple
%   functions with ORDER specifying tensor layout ('first' or 'last').
%
%   This function serves as an accessor for tensor-based algorithms like
%   DMRG_CROSS. It computes values of univariate function(s) at specific
%   indices of the (joint) QTT tensor defined on a uniform grid.
%
%   See also: POLY2QTT, DMRG_CROSS

    if nargin < 6
        M = 1;
    end
    if nargin < 5
        order = '';
    end

    n = 2^d;
    h = (dom(2) - dom(1)) / n;

    if M == 1
        % --- Single Function Mode ---
        grid_multi_indices = ind;
    else
        % --- Multi-Function (Vectorized) Mode ---
        if strcmpi(order, 'first')
            % Tensor is M x 2 x ... x 2. ind is [j, i1, ..., id]
            poly_indices = ind(:, 1);
            grid_multi_indices = ind(:, 2:d+1);
        elseif strcmpi(order, 'last')
            % Tensor is 2 x ... x 2 x M. ind is [i1, ..., id, j]
            poly_indices = ind(:, d+1);
            grid_multi_indices = ind(:, 1:d);
        else
            error('Invalid order specified. Use ''first'' or ''last''.');
        end
    end

    % Convert base-2 multi-index to linear index 'k' (little-endian)
    linear_indices = 1 + sum((grid_multi_indices - 1) .* 2.^(0:d-1), 2);

    % Calculate grid points 's_k' from linear indices.
    grid_points = dom(1) + (linear_indices - 0.5) * h;

    % Evaluate function(s) at the grid points.
    all_vals_at_points = fun(grid_points); % size: length(grid_points) x M

    if M == 1
        fun_tensor = all_vals_at_points;
    else
        % Select the values for the specific functions requested by poly_indices.
        num_requested_points = size(ind, 1);
        selection_indices = sub2ind(size(all_vals_at_points), (1:num_requested_points)', poly_indices);
        fun_tensor = all_vals_at_points(selection_indices);
    end
end
