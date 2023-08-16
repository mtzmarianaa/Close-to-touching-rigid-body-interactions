function mat = evaluateRInterpolated(xDist, listPrecomputedR)
% *evaluateRInterpolated* builds the R matrix interpolated from precomputed R matrices 
% evaluated at interpolation nodes.
%
% Syntax: mat = evaluateRInterpolated(xDist, listPrecomputedR)
%
% Input:
%   xDist - distance between two discs
%   listPrecomputedR - list of precomputed R matrices for different
%                                 distances
%
% Output:
%   mat - interpolated R matrix
%
% To do: change inteprolation nodes to log Chebyshev.
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

% Given a k matrices of sizes nxn (matrices evaluated at different
% distances
% We interpolate using those matrices to evaluate K at xDist

n = size(listPrecomputedR, 1);
mat = zeros( n,n ); % Initialize matrix
k = size(listPrecomputedR, 3);

% Compute weights for barycentric Lagrange interpolation
w = lege.barywts(k); % kx1 vector
x = lege.exps(k); % Lagrange nodes

W = w./(xDist - x);
denom = sum(W);

A = reshape( listPrecomputedR, [], k);

mat = A*W;
mat = reshape( mat, n, []);
mat = mat./denom;

end