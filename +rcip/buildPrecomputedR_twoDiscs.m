function [listPrecomputedR] = buildPrecomputedR_twoDiscs(geom0, pClose0, listK22_inv, nRefDist, kDist)
% *buildPrecomputedR_twoDiscs* builds the precomputed R matrix for
% different distances and saves them into a cell array. 
%
% Syntax: buildPrecomputedR_twoDiscs(geom0, pClose0, listK22_inv)
%              buildPrecomputedR_twoDiscs(geom0, pClose0, listK22_inv, nRefDist)
%              buildPrecomputedR_twoDiscs(geom0, pClose0, listK22_inv, nRefDist, kDist)
%
% Input:
%   geom0 - initial geometry of the discs for distance0 (for the
%                interpolation on the distance between two discs)
%   pClose0 - information of the close-to-toching region of the discs with
%                   initial distance 
%   listK22_inv - list of inverses for the close-to-touching region,
%                      $K_{22}^{-1}$ for different distances.
%
% Optional input:
%   nRefDist - Number of chunks to use in the discretization for the distance
%   kDist - number of discretization points to use for the distance
%
% Output:
%   listPrecomputedR - list of precomputed R matrices for desired
%                                 distances.
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)


if(nargin<4)
    nRefDist = 28;
end

if(nargin<5)
    kDist = 16;
end

ds = discs(geom0, pClose0);

% Build the breakpoints for the distance, compute the chaning center for
% disc 2 (moving closer to disc 1)
ks = 0:nRefDist;
breakPoints_dist = geom0.Rs(1) + 1./2.^ks*( norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2) ) + geom0.Rs(2);
cHat2 = ( geom0.ctrs(:,2) - geom0.ctrs(:,1))./( norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) )*breakPoints_dist + geom0.ctrs(:, 1);

% For the Legendre nodes for the distances
xDist= lege.exps(kDist);

% Initialize the cell array and the cell item
listPrecomputedR = cell(1, nRefDist);
cellItem = zeros(ds.listCoarseGammas(1).npt, ds.listCoarseGammas(1).npt, kDist);
cInter = zeros(2, kDist);

geom.Rs = geom0.Rs;
geom.nBreakPoints = geom0.nBreakPoints;

% Fill the cell array
for i=1:nRefDist
    % Legendre nodes center of disc 2
    a = cHat2(1,i);
    b = cHat2(1, i+1);
    cInter(1, :) = a + (b-a)*(xDist + 1)/2;
    a = cHat2(2,i);
    b = cHat2(2, i+1);
    cInter(2, :) = a + (b-a)*(xDist + 1)/2;
    for j=1:kDist
        K22_inv = listK22_inv{(i-1)*kDist + j};
        % Compute R for each of these nodes
        geom.ctrs(:, 1) = geom0.ctrs(:, 1);
        geom.ctrs(:, 2) = cInter(:, j); % New center
        ds = discs(geom, pClose0);
        nRef = floor(ds.listGammas(1).nch/4 - 2);
        P = rcip.prol_dyadic(ds.listCoarseGammas(1).k, nRef);
        P = blkdiag(P, P);
        cellItem(:, :, j) = rcip.buildR(ds.listCoarseGammas(1), ds.listGammas(1), K22_inv, P);
    end
    listPrecomputedR{i} = cellItem;
end


end