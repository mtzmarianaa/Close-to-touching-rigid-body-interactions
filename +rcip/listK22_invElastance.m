function listK22_inv = listK22_invElastance(geom0, pClose0, nRefDist, kDist)
% Computes the list of K22 inverses of two discs for the elastance
% problem (as the centers move closer and we dyadically refine on the
% distance)

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

kern = @(s,t) krns.DL_kern(s,t);
flagFunction = @(i, ds) flagFunctionGamma1(i, ds);

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
        nB = [0, ds.listGammas(1).npt];
        M = dsc.buildMelastance(ds, ds.listGammas, nB, flagFunction);
        matOffSet = 0.5*eye(ds.listGammas(1).npt) + M;
        K22 = chunkermat(ds.listGammas(1), kern, opts2) + matOffSet;
        listK22_inv{(i-1)*kDist + j} = inv(K22);
%         fprintf("%3.0e :  %8.0e  x  %8.0e \n", (i-1)*kDist + j, size(K22));
%         (i-1)*kDist + j
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