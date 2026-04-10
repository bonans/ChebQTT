function subs = lin2sub_all2(d, lin_inds)
%LIN2SUB_ALL2 Convert linear indices to binary multi-indices.
%   SUBS = LIN2SUB_ALL2(D, LIN_INDS) converts linear indices to
%   D-dimensional binary multi-indices.
%
%   LIN_INDS : Array of linear indices (1 to 2^D)
%   SUBS     : Array of size numel(LIN_INDS) x D with values in {1,2}
%
%   See also: TUCKQTT_EVAL

    if isvector(lin_inds)
        N = length(lin_inds);
        sz = N;
    else
        sz = size(lin_inds);
        N  = numel(lin_inds);
    end

    % zero-based linear indices
    lin0 = lin_inds(:) - 1;

    % allocate output
    subs = zeros([N d]);

    % extract bits
    for j = 1:d
        subs(:, j) = bitget(lin0, j);
    end
    subs = subs + 1;

    % reshape back to tensor form
    subs = reshape(subs, [sz d]);
end
