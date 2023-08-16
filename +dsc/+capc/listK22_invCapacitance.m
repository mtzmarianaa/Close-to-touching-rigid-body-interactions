function listK22_inv = listK22_invCapacitance(geom0, pClose0, nRefDist, kDist)
% *listK22_invCapacitance* computes the list of $K_{22}^{-1}$ of two discs
% for the capacitance problem. As discs move closer and closer together.
%
% Syntax: listK22_inv = listK22_invCapacitance(geom0, pClose0)
%              listK22_inv = listK22_invCapacitance(geom0, pClose0, nRefDist)
%              listK22_inv = listK22_invCapacitance(geom0, pClose0, nRefDist, kDist)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   geom0 - initial geometry of the discs for distance0 (for the
%                interpolation on the distance between two discs)
%   pClose0 - information of the close-to-toching region of the discs with
%                   initial distance 
%
% Optional input:
%   nRefDist - Number of chunks to use in the discretization for the distance
%   kDist - number of discretization points to use for the distance
%
% Output:
%   listK22_inv - list of inverses for the close-to-touching region,
%                      $K_{22}^{-1}$ for different distances.
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
xDist = lege.exps(kDist);

% Initialize the cell array and the cell item
listK22_inv = cell(1, kDist*nRefDist);

geom.Rs = geom0.Rs;
geom.nBreakPoints = geom0.nBreakPoints;

% Options for on boundary evaluations
opts2 = [];
opts2.adaptive_correction = true;

kern = @(s,t) krns.DLplusSL(s,t);

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
        geom.ctrs(:, 1) = geom0.ctrs(:, 1);
        geom.ctrs(:, 2) = cInter(:, j); % New center
        ds = discs(geom, pClose0);
        matOffSet = 0.5*eye(ds.listGammas(1).npt);
        K22 = chunkermat(ds.listGammas(1), kern, opts2) + matOffSet;
        listK22_inv{(i-1)*kDist + j} = inv(K22);
%         fprintf("%3.0e :  %8.0e  x  %8.0e \n", (i-1)*kDist + j, size(K22));
%         (i-1)*kDist + j
    end
end

end