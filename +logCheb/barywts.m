function w = barywts(x)
% *barywts* returns scaled barycentric weights for points x (NOT WEIGHTED)
% Makes use of chebfun! q can be either the number of points or the points
% themselves. Similar to *CHEBFUN*!
%
% Syntax: w = logCheb.barywts(x)
% Input:
%   x - nodes
%
% Output:
%   w - Barycentric weights
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

k = length(x);

xx = bsxfun(@minus,x,x.');
xx(1:k+1:end) = 1;

w = prod(xx,1);
w = w(1)./w(:);


end