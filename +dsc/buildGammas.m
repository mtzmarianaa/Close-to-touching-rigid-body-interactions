function [listGammas, neisMapClose] = buildGammas(geom, pClose, I, nDiscs)
% Builds the list of gammas, the close to touching part. Also build
% neisMapClose, the map for the neighbors. Not user friendly

pref = chunkerpref();
cparams = [];
cparams.ta = -2*pi;
cparams.tb = 4*pi;

nGammas = size(I, 1); % Number of close to touching parts
listGammas(1, nGammas) = chunker(); % Initialize the list of chunker objects
neisMapClose = []; % Initialize struct

nchold = 0; % Initialize off set
neisClose = [];
nCloseT = 0;

% Initialize the neisClose object
for i=1:nDiscs
    neisClose(i).map = -1*ones(2*pClose(i).nClose, 2 );
    neisClose(i).currentInd = 1;
    nCloseT = nCloseT + 2*pClose(i).nClose;
end

% Start building
for r=1:nGammas

    % For Di
    i = I(r, 1);
    currentInd_i = neisClose(i).currentInd
    thetaReg_i = pClose(i).thetasReg(currentInd_i);
    R_i = geom.Rs(i);
    ctr_i = geom.ctrs(:, i);
    theta_ji = I(r, 2);
    % For Ds
    s = I(r, 3);
    currentInd_s = neisClose(s).currentInd;
    thetaReg_s = pClose(s).thetasReg(currentInd_s);
    R_s = geom.Rs(s);
    ctr_s = geom.ctrs(:, s);
    theta_ks = I(r, 4);

    % Compute distance, alpha_0
    dist = dist_is(ctr_i, R_i, theta_ji, ctr_s, R_s, theta_ks);
    
    % Start building the part of gamma_r from Di, add information to
    % neisClose(i)
    alpha_0i = minAlpha(thetaReg_i, R_i, dist);
    plusReg_i = fanReg(alpha_0i, theta_ji); % breakpoints
    d_i = @(t) disc(t, center = ctr_i, radius = R_i);
    chnkr_gamma_r1 = chunkerfuncbreakpoints( d_i, plusReg_i, cparams, pref, false );
    plusReg_i = mod(plusReg_i, 2*pi);
    startPoint_i = plusReg_i(1, 1);
    nStart_i = 1 + nchold;
    nchold = nchold + chnkr_gamma_r1.nch;
    endPoints_i = plusReg_i(1, end);
    nEnd_i = nchold;
    neisClose(i).map(currentInd_i, 1) = startPoint_i;
    neisClose(i).map(currentInd_i, 2) = nStart_i;
    neisClose(i).map(currentInd_i + 1, 1) = endPoints_i ;
    neisClose(i).map(currentInd_i + 1, 2) = nEnd_i;
    neisClose(i).currentInd = neisClose(i).currentInd + 1;


    % Start building the part of gamma_r from Ds, add information to
    % neisClose(s)
    alpha_0s = minAlpha(thetaReg_s, R_s, dist);
    plusReg_s = fanReg(alpha_0s, theta_ks); % breakpoints
    d_s = @(t) disc(t, center = ctr_s, radius = R_s);
    chnkr_gamma_r2 = chunkerfuncbreakpoints( d_s, plusReg_s, cparams, pref, false );
    plusReg_s = mod(plusReg_s, 2*pi);
    startPoint_s = plusReg_s(1, 1);
    nStart_s = 1 + nchold;
    nchold = nchold + chnkr_gamma_r2.nch;
    endPoints_s = plusReg_s(1, end);
    nEnd_s = nchold;
    neisClose(s).map(currentInd_s, 1) = startPoint_s;
    neisClose(s).map(currentInd_s, 2) = nStart_s;
    neisClose(s).map(currentInd_s + 1, 1) = endPoints_s ;
    neisClose(s).map(currentInd_s + 1, 2) = nEnd_s;
    neisClose(s).currentInd = neisClose(s).currentInd + 1;

    % Merge the two parts of gamma_r, add them to the list
    listGammas(1, r) = merge( [chnkr_gamma_r1, chnkr_gamma_r2] );

end

% Now sort neisClose(:).map, initialize neisMapClose
neisMapClose = zeros(nCloseT, 1);
rowLims = zeros(nDiscs + 1, 1);

for i=1:nDiscs
    neisClose(i).map = sortrows(neisClose(i).map);
    rowLims(i+1, 1) = rowLims(i, 1) + 2*pClose(i).nClose;
end

% Fill the matrix neisMapClose

for i=1:nDiscs
    neisMapClose(rowLims(i, 1)+1: rowLims(i+1, 1), 1) = neisClose(i).map(:, 2);
end


end


function dist = dist_is(ctr_i, R_i, theta_ji, ctr_s, R_s, theta_ks)
% Computes the minimum distance between disc i and disc s
% (distance between Di(theta_ji) and Ds(theta_ks)

x_ji = disc(theta_ji, center = ctr_i, radius = R_i);
x_ks = disc(theta_ks, center = ctr_s, radius = R_s);

dist = norm(x_ji - x_ks);

end

function alpha_0 = minAlpha(thetaReg, R, dist)
% Computes minimum refinement needed

alpha_0 = thetaReg;

while( R*alpha_0 > dist )
    alpha_0 = 0.5*alpha_0;
end


end



function plusReg = fanReg(alpha_0, theta)

L = log(pi/(6*alpha_0))/log(2);
L = max(0, round(L)); % just in case

plusReg = zeros(1, L+1);

for i=0:L
    plusReg(1, (i+1)) = alpha_0*2^i;
end

negPlusReg = -flip(plusReg);

plusReg = [negPlusReg + theta theta plusReg + theta];
%plusReg = mod(plusReg, 2*pi); % mod 2pi for sorting purposes

end
