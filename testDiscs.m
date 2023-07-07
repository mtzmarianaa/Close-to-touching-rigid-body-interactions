% Test discs


% 3 discs

geom = [];
geom.ctrs = [0 1.75 0; 0 0 1.75];
geom.Rs = [0.75; 0.75; 0.75];
geom.nBreakPoints = [20; 20; 20];
pClose = [];
pClose(1).data = [0 2 1; pi/2 3 1];
pClose(1).nClose = 2;
pClose(1).thetasReg = [pi/6 pi/6];
pClose(2).data = [pi, 1, 1];
pClose(2).nClose = 1;
pClose(2).thetasReg = pi/6;
pClose(3).data = [3*pi/2 1 2];
pClose(3).nClose = 1;
pClose(3).thetasReg = pi/6;

nDiscs = 3;

% Build the necessary things


[I, nGammas, pClose] = dsc.buildI(pClose, nDiscs);
[listGammas, neisMapClose] = dsc.buildGammas(geom, pClose, I, nDiscs);
chnkrsGammas = merge(listGammas);
[gamma0, neisMapFar, listFarChunks] = dsc.buildGamma0(geom, pClose, I, nDiscs);

% Add the missing neis
indMissingClose = find(~chnkrsGammas.adj);
indMissingFar = find(~gamma0.adj);
% Merge
chnkrs = merge([gamma0, listGammas]);
chnkrs .adj(2*gamma0.nch + indMissingClose) = neisMapFar;
chnkrs.adj(indMissingFar) = gamma0.nch + neisMapClose;

indGamma = gamma0.nch + 1;



% Plot things
figure()
plot(gamma0, '-o')
axis equal
title("Far section of the curve")


figure()
plot(chnkrsGammas, '-*')
axis equal
title("Near section of the curve")



figure()
plot(chnkrs, '-.')
axis equal
title("Far and close sections of the curved combined")

