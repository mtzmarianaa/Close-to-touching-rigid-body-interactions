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
k = size(listPrecomputedR, 3);

% Compute weights for barycentric Lagrange interpolation or log Chebyshev
if(typeNodes == 'l')
    w = lege.barywts(k); % kx1 vector
    x = lege.exps(k); % Lagrange nodes
else
    [x, w] = logCheb.exps(k);
end

% If we have hit a node
hitNode = abs(xDist - x)<1e-13;
if any(hitNode)
    ind = find(hitNode);
    mat = listPrecomputedR(:, :, ind);
    return
end

W = w./(xDist - x);
denom = sum(W);

A = reshape( listPrecomputedR, [], k);

mat = A*W;
mat = reshape( mat, n, []);
mat = mat./denom;


% For debugging, something is wrong the interpolation for log Cheby

u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);
uk = {u1, u2};
kern = kernel('lap', 'c', [1.0, 1.0]);
opts2 = [];
opts2.adaptive_correction = true;

d = 0.1 + (0.1 - 1e-12)*(xDist + 1)/2;

pClose = [];
pClose(1).data = [0 2 1];
pClose(1).nClose = 1;
pClose(1).thetasReg = pi/6;
pClose(2).data = [pi, 1, 1];
pClose(2).nClose =1;
pClose(2).thetasReg = pi/6;

ctrs = [0 1.5 ;0 1.5+d]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
geom.ctrs = ctrs;

ds = discs(geom, pClose);

nRef = floor(ds.listGammas(1).nch/4 - 2);
P = rcip.prol_dyadic(ds.listCoarseGammas(1).k, nRef);
P = blkdiag(P, P);
matOffSet = 0.5*eye(ds.listGammas(1).npt);
K22 = chunkermat(ds.listGammas(1), kern, opts2) + matOffSet;
K22_inv = inv(K22);
trueMat = rcip.buildR(ds.listCoarseGammas(1), ds.listGammas(1), K22_inv, P);


n = size(listPrecomputedR, 1);
k = size(listPrecomputedR, 3);
xDistVec = linspace(-1, 1, 1000)';
r = 25;
c = 25;
interp = zeros(100, size(r, 1));
trueSubMat = zeros(100, size(r, 1));
for i=1:1000
    W = w./(xDistVec(i, 1) - x);
    denom = sum(W);
    A = reshape( listPrecomputedR, [], k);
    mat = A*W;
    mat = reshape( mat, n, []);
    mat = mat./denom;
    interp(i, :) = mat(r, c);
    % For the true solution
    d = 0.1 + (-0.1 + 1e-12)*(xDistVec(i, 1) + 1)/2;
    ctrs = [0 1.5 ;0 1.5+d];
    geom.ctrs = ctrs;
    ds = discs(geom, pClose);
    nRef = floor(ds.listGammas(1).nch/4 - 2);
    P = rcip.prol_dyadic(ds.listCoarseGammas(1).k, nRef);
    P = blkdiag(P, P);
    matOffSet = 0.5*eye(ds.listGammas(1).npt);
    K22 = chunkermat(ds.listGammas(1), kern, opts2) + matOffSet;
    K22_inv = inv(K22);
    trueMat = rcip.buildR(ds.listCoarseGammas(1), ds.listGammas(1), K22_inv, P);
    trueSubMat(i, :) = trueMat(r,c);
end


end