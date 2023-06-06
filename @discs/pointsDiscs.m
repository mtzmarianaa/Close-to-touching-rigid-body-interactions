function [r, d, d2, endDiscs] = pointsDiscs(t,varargin)

% Get points on discs
% return position, first, second derivatives of parameterized discs,
% and the indices of the last point on each disc
%
% Inputs:
% t - points (in [0,2pi] to evaluate these quantities
%     (note that since we are dealing with identical discs, if we have 
%      nk discs in total then (j-1)*2pi/nk <= t <= j*2pi/nk for the
%      j-th disc)
%
% Optional inputs:
% R - float, radius of the discs (1)
% ctrs - float(2, ndiscs), x,y coordinates of center of discs ( [0,0] )

ndiscs = 1;
R = 1.0;
ctrs = [0;0];

if nargin > 1 && ~isempty(varargin{1})
    R = varargin{2};
end
if nargin > 1 && ~isempty(varargin{2})
    ctrs = varargin{2};
    ndiscs = size(ctrs, 2);
end

% Make sure we have the correct amount of points 
nPoints = length(t);

if(mod(nPoints, ndiscs) ~= 0)
    nPoints = ceil(nPoints/ndiscs)*ndiscs;
    t = linspace(0, 2*pi, nPoints);
end

endDiscs = zeros(1, ndiscs);
r = zeros(2, nPoints);
d = zeros(2, nPoints);
d2 = zeros(2, nPoints);

% Compute r, d, d2

ct = cos(ndiscs*t);
st = sin(ndiscs*t);

% fill in by hand the parametrization for the first circle

indIncluded = find( t<= 2*pi/ndiscs );
endDiscs(1) = indIncluded(end);
r(:, 1:endDiscs(1) ) = [R*ct(1:endDiscs(1)) + ctrs(1, 1); R*st(1:endDiscs(1)) + ctrs(2, 1)];
d(:, 1:endDiscs(1) ) = [-st(1:endDiscs(1)); ct(1:endDiscs(1))];
d2(:, 1:endDiscs(1) ) = [-ct(1:endDiscs(1)); -st(1:endDiscs(1))];
 
for j=2:ndiscs
    indIncluded = find(t> (j-1)*2*pi/ndiscs & t<=j*2*pi/ndiscs);
    endDiscs(j) = indIncluded(end); % end of disc j
    r(:, endDiscs(j-1):endDiscs(j) ) = [R*ct(endDiscs(j-1):endDiscs(j)) + ctrs(1, j); R*st(endDiscs(j-1):endDiscs(j)) + ctrs(2, j)];
    d(:, endDiscs(j-1):endDiscs(j) ) = [-st(endDiscs(j-1):endDiscs(j)); ct(endDiscs(j-1):endDiscs(j)) ];
    d2(:, endDiscs(j-1):endDiscs(j) ) = [-ct(endDiscs(j-1):endDiscs(j)); -st(endDiscs(j-1):endDiscs(j)) ];
end

end