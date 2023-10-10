function mat = evaluateRInterpolated(xDist, listPrecomputedR, typeNodes)
% *evaluateRInterpolated* builds the R matrix interpolated from precomputed R matrices 
% evaluated at interpolation nodes.
%
% Syntax: mat = evaluateRInterpolated(xDist, listPrecomputedR)
%              mat = evaluateRInterpolated(xDist, listPrecomputedR, typeNodes)
%
% Input:
%   xDist - distance between two discs
%   listPrecomputedR - list of precomputed R matrices for different
%                                 distances
%
% Optional input:
%   typeNodes - 'l' for Legendre nodes, 'logc' for log Chebyshev (auto 'l')
%
% Output:
%   mat - interpolated R matrix
%
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

if(nargin < 3)
    typeNodes = 'logc';
else
    typeNodes = lower(typeNodes);
end

% Given a k matrices of sizes nxn (matrices evaluated at different
% distances
% We interpolate using those matrices to evaluate K at xDist

n = size(listPrecomputedR, 1);
mat = zeros( n,n ); % Initialize matrix
k = size(listPrecomputedR, 3);

% Compute weights for barycentric Lagrange interpolation
if(typeNodes == 'l')
    w = lege.barywts(k); % kx1 vector
    x = lege.exps(k); % Lagrange nodes
else
    [x, w] = exps(k);
end

W = w./(xDist - x);
denom = sum(W);

A = reshape( listPrecomputedR, [], k);

mat = A*W;
mat = reshape( mat, n, []);
mat = mat./denom;

end