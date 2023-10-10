function w = barywts(q, a, b)
% *barywts* returns scaled barycentric weights for the logChebyshev points
% Makes use of chebfun! q can be either the number of points or the points
% themselves. Makes use of *CHEBFUN*!
%
% Syntax: w = logCheb.barywts(q)
%              w = logCheb.barywts(q, a, b)
% Input:
%   q - either number of nodes or nodes themselves
%
% Optional input:
%   a - start of the interval if q is number of nodes (auto -1)
%   b - end of interval if q is number of nodes (auto 1)
%
% Output:
%   w - Barycentric weights
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

if(nargin < 2)
    a = -1;
    b = 1;
end

if length(q) == 1
    % Means that we just have k, the number of points
    [~, w] = logCheb.exps(q, a, b);
else
    % q is a vector
    w = baryWeights(q);
end


end