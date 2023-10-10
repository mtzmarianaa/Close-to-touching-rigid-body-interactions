function [x, w] = exps(k, a, b)
% logCheb.exps gets the nodes and the weights for log(Chebyshev)
% Makes use of *CHEBFUN*!
%
% Syntax: x = logCheb.barywts(k)
%              x = logCheb.barywts(k, a, b)
%              [x, w] = logCheb.barywts(k)
%              [x, w] = logCheb.barywts(k, a, b)
% Input:
%   k - either number of nodes
%
% Optional input:
%   a - start of the interval if q is number of nodes (auto -1)
%   b - end of interval if q is number of nodes (auto 1)
%
% Output:
%   x - log Chebyshev nodes
%   w - Barycentric weights
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

if(nargin<2)
    a = -1;
    b = 1;
end

% Generate the log(Cheby) points from a to b
x = chebpts(k+1); % Chebyshev points of the second kind
x = (x + 1)/2;
x = x(2:end, :);
x = log(x);
x = a + (b-a)*(x - x(1))/(x(end) - x(1));

if(nargout>1)
    w = logCheb.barywts(x); % Barycentric weights from chebfun
end


end