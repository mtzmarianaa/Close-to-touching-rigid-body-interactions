function [gamma0, neisMapFar, listFarChunks] = buildGamma0(geom, pClose, I, nDiscs)
% Function that builds the far intervals (for gamma 0). Also builds
% neisMapFar, the map for the neighbors. Not user friendly. Same idea as
% buildGammas

pref = chunkerpref();
cparams = [];
cparams.ta = -2*pi;
cparams.tb = 4*pi;

nGammas = size(I, 1); % Number of close to touching parts
neisMapFar = zeros(4*nGammas, 4); % With this we are going to sort the neighbor indices
farInfo = [];
offSetNeisMap = zeros(nDiscs + 1, 1);
neisMapFar(:, 1) = 1:4*nGammas;
listFarChunks(1, nDiscs) = chunker();

% Initialize farInfo
for i=1:nDiscs
    farInfo(i).breakpointsClose = zeros(pClose(i).nClose, 3);
    farInfo(i).farIntervals = zeros(pClose(i).nClose, 2);
    farInfo(i).lenFarIntervals = zeros(pClose(i).nClose, 1);
    farInfo(i).currentIndex = 1;
    offSetNeisMap(i + 1) = offSetNeisMap(i) + 2*pClose(i).nClose;
end

% Build breakpointsClose and start filling in neisMapFar
for r=1:nGammas
    r4 = 4*r - 3;
    % For Di
    i = I(r, 1);
    currentInd_i = farInfo(i).currentIndex;
    thetaReg_i = pClose(i).thetasReg(currentInd_i);
    theta_ji = I(r, 2);
    br_i = [ theta_ji - thetaReg_i; theta_ji;  theta_ji + thetaReg_i];
    farInfo(i).breakpointsClose(currentInd_i, :) = br_i;
    farInfo(i).currentIndex = currentInd_i  + 1;
    % For Ds
    s = I(r, 3);
    currentInd_s = farInfo(s).currentIndex;
    thetaReg_s = pClose(s).thetasReg(currentInd_s);
    theta_ks = I(r, 4);
    br_s = [ theta_ks - thetaReg_s; theta_ks;  theta_ks + thetaReg_s];
    farInfo(s).breakpointsClose(currentInd_s, :) = br_s;
    farInfo(s).currentIndex = currentInd_s  + 1;
    % Add this information to neisMapClose
    neisMapFar(r4, 2) = br_i(1);
    neisMapFar(r4+1, 2) = br_i(3);
    neisMapFar(r4+2, 2) = br_s(1);
    neisMapFar(r4+3, 2) = br_s(3);
    neisMapFar(r4, 3) = i;
    neisMapFar(r4+1, 3) = i;
    neisMapFar(r4+2, 3) = s;
    neisMapFar(r4+3, 3) = s;
end

% sort
%neisMapFar(:, 2) = mod(neisMapFar(:, 2), 2*pi);
neisMapFar = sortrows(neisMapFar, [3, 2]);

% Build the far pieces
listFarPieces(1, nDiscs) = chunker();

% Build far intervals, compute their lengths
for r=1:nDiscs
    % Compute ends in neisMapFar
    startN = offSetNeisMap(r) + 1;
    endN = offSetNeisMap(r + 1);
    % Sort farInfo
    farInfo(r).breakpointsClose = sort(farInfo(r).breakpointsClose, 2); % Order per closest angle
    % Switch first and last row in this section of offSetNeisMap
    fR = neisMapFar(startN, :);
    neisMapFar(startN:(endN - 1), :) = neisMapFar( (startN + 1):endN, : );
    neisMapFar(endN, :) = fR;
    % Build far intervals and lengths
    for j=1:(pClose(r).nClose -1)
        a = farInfo(r).breakpointsClose(j, 3);
        b = farInfo(r).breakpointsClose(j+1, 1);
        if( b < a )
            b = b + 2*pi;
        end
        farInfo(r).farIntervals(j, :) = [ a b ];
        farInfo(r).lenFarIntervals(j) = b - a;
    end
    % For the last one
    a = farInfo(r).breakpointsClose(end, 3);
    b = farInfo(r).breakpointsClose(1, 1);
    if( b < a)
        b = b + 2*pi;
    end
    farInfo(r).farIntervals(end, :) = [a b];
    farInfo(r).lenFarIntervals(end) = b - a;
end


% Compute how many breakpoints we need in each of the far intervals based
% on nBreakpoints and the lengths of the intervals, build a chunker object
% per disc

currentNei = 0;

for r=1:nDiscs
    startN = offSetNeisMap(r);
    d_r = @(t) disc(t, center = geom.ctrs(:, r), radius = geom.Rs(r));
    % Compute how many breakpoints we have
    nBreakpoints_r = geom.nBreakPoints(r) - 3*pClose(i).nClose;
    lenFarSection = sum( farInfo(r).lenFarIntervals ); % Total length of the far section
    for j=1:pClose(r).nClose
        % Compute how many breakpoints to put in this interval
        nB = round( nBreakpoints_r*farInfo(r).lenFarIntervals(j)/lenFarSection );
        a = farInfo(r).farIntervals(j, 1);
        b = farInfo(r).farIntervals(j, 2);
        ts = linspace( a, b, max(2, nB) );
        % Build this piece of the chunker
        chnkr_rj = chunkerfuncbreakpoints(d_r, ts, cparams, pref, false);
        if( j == 1)
            listFarChunks(1, r) = chnkr_rj;
        else
            listFarChunks(1, r) = merge( [listFarChunks(1, r), chnkr_rj] );
        end
        % Add the neighbor index
        neisMapFar(startN + 2*j - 1, 4) = currentNei + 1;
        currentNei = currentNei + chnkr_rj.nch;
        neisMapFar(startN + 2*j, 4) = currentNei;
    end
    % Flip
%     lR = neisMapFar( offSetNeisMap(r+1), 4);
%     neisMapFar(startN+2:(offSetNeisMap(r+1) ), 4  ) = neisMapFar(startN+1:(offSetNeisMap(r+1) - 1), 4  );
%     neisMapFar(startN + 1, 4) = lR;
end


% Combine all of them

gamma0 = merge( listFarChunks );

neisMapFar = sortrows(neisMapFar);
neisMapFar = neisMapFar(:, 4);
end