% Test discs

addpaths_loc();
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



% Directly with the builder

ds = discs(geom, pClose);
