function points = I2unitI(points, I)
%I2UNITI Transform points from interval I to the unit interval [-1,1].
%   POINTS = I2UNITI(POINTS, I) maps points from the interval I = [a,b]
%   to the standard interval [-1,1].
%
%   See also: CHEBPOLYEVAL

points = 2 * (points - I(1)) / (I(2) - I(1)) - 1;

% Clamp the results to be strictly within [-1, 1]
points = max(-1, min(1, points));
end
