function listK22_inv = listK22_invElastance(geom0, pClose0, typeNodes, nRefDist, kDist)
% *listK22_invCapacitance* computes the list of $K_{22}^{-1}$ of two discs
% for the elastance problem. As discs move closer and closer together.
%
% Syntax: listK22_inv = listK22_invElastance(geom0, pClose0)
%              listK22_inv = listK22_invElastance(geom0, pClose0, nRefDist)
%              listK22_inv = listK22_invElastance(geom0, pClose0, nRefDist, kDist)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   geom0 - initial geometry of the discs for distance0 (for the
%                interpolation on the distance between two discs)
%   pClose0 - information of the close-to-toching region of the discs with
%                   initial distance 
%   typeNodes - 'l' for Legendre nodes, 'logc' for log Chebyshev
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

if(isequal(typeNodes, 'logc'))
    kDist = 16;
    nRefDist = 1;
end

% Build the breakpoints for the distance, compute the changing center for
% disc 2 (moving closer to disc 1)
ks = 0:nRefDist;
if(isequal(typeNodes, 'l'))
    breakPoints_dist = geom0.Rs(1) + 1./2.^ks*( norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2) ) + geom0.Rs(2);
else
    breakPoints_dist = geom0.Rs(1) + geom0.Rs(2) + [0.1 1e-12];
end

cHat2 = ( geom0.ctrs(:,2) - geom0.ctrs(:,1))./( norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) )*breakPoints_dist + geom0.ctrs(:, 1);

% For the Legendre nodes for the distances
if(isequal(typeNodes, 'l'))
    xDist = lege.exps(kDist);
else
    % For log Chebyshev of the distances
    xDist = logCheb.exps(kDist, breakPoints_dist(1), breakPoints_dist(2));
end

% Initialize the cell array and the cell item
listK22_inv = cell(1, kDist*nRefDist);

geom.Rs = geom0.Rs;
geom.nBreakPoints = geom0.nBreakPoints;

% Options for on boundary evaluations
opts2 = [];
opts2.adaptive_correction = true;

kern = kernel('lap', 'd');
flagFunction = @(i, ds) flagFunctionGamma1(i, ds);


if(isequal(typeNodes, 'l'))
    % Fill the cell array for Legendre nodes
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
            nB = [0, ds.listGammas(1).npt];
            M = dsc.elst.buildMelastance(ds, ds.listGammas, nB, flagFunction);
            matOffSet = 0.5*eye(ds.listGammas(1).npt) + M;
            K22 = chunkermat(ds.listGammas(1), kern, opts2) + matOffSet;
            listK22_inv{(i-1)*kDist + j} = inv(K22);
        end
    end
else
    % Fill the cell array for log Chebyshev
    geom.ctrs(:, 1) = geom0.ctrs(:, 1);
    for j=1:kDist
        geom.ctrs(1, 2) = xDist(j); % New center for second disc
        ds = discs(geom, pClose0);
        nB = [0, ds.listGammas(1).npt];
        M = dsc.elst.buildMelastance(ds, ds.listGammas, nB, flagFunction);
        matOffSet = 0.5*eye(ds.listGammas(1).npt) + M;
        K22 = chunkermat(ds.listGammas(1), kern, opts2) + matOffSet;
        listK22_inv{j} = inv(K22);
    end
end
end


function flagGamma1 = flagFunctionGamma1(i, ds)

flagGamma1 = zeros(ds.listGammas(1).nch, 1);

% Look for the close to touching parts using I and add
% those who belong to disc i
disc_i = ds.I(1, 1);
disc_s = ds.I(1, 3);
k = 1;
while( (disc_i <= i || disc_s <= i) && k <= ds.nGammas )
    disc_i = ds.I(k, 1);
    disc_s = ds.I(k, 3);
    if(disc_i == i)
        % Meaning the first part of this gamma is on disc i
        flagGamma1( (ds.listGammas(1).nch/2 + 1):end, 1 ) = 1;
    end
    if( disc_s == i)
        % Meaning that the second part of this gamma is on disc i
        flagGamma1( 1:(ds.listGammas(1).nch/2), 1 ) = 1;
    end
    % Update
    k = k+1;
end

end