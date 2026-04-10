function r = TTrankMean(tt)
%TTRANKMEAN Compute the geometric mean of TT ranks.
%   R = TTRANKMEAN(TT) computes the geometric mean of the TT ranks
%   of the tensor TT.
%
%   See also: TT_TENSOR

r = rank(tt);
r = sqrt(mean(r(1:end-1).*r(2:end)));
end
