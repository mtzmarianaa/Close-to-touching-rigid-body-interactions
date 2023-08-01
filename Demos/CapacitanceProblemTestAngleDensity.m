%%% Capacitance problem test solution density vs angle

addpaths_loc();
clear all
close all
clc

u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);

uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 1.500005 ;0 0]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [20; 20];

%%%%%%%%%% UNIFORM PANELS
% Define points on surface
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;

pClose = [];
pClose(1).data = [0 2 1];
pClose(1).nClose = 1;
pClose(1).thetasReg = pi/6;
pClose(2).data = [pi, 1, 1];
pClose(2).nClose =1;
pClose(2).thetasReg = pi/6;

xcoordCtr2 = 1.5 + [ 1e-1 1e-3 1e-5];
nTest = length(xcoordCtr2);

figure()
hold on

colors = [0 23 255; 111 0 203; 60 228 225];
colors = colors./255;

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    geom.ctrs = ctrs;
    ds2 = discs(geom, pClose);
    [q, sigma] = capacitanceProblem(ds2, uk);
    % Get the angles
    flag2 = dsc.flagnDisc(2, ds2);
    flag2_points = repmat(flag2, 1, ds2.chnkrs.k);
    flag2_points = flag2_points';
    flag2_points = logical( flag2_points(:) );
    points2 = reshape( ds2.chnkrs.r, 2, [] );
    points2 = points2(:, flag2_points);
    sigma2 = sigma(flag2_points);
    addPiTh = points2(2, :) < 0;
    thetas2 = acos(  (points2(1, :) - ctrs(1, 2) )./( Rs(2) ) );
    thetas2(addPiTh) = 2*pi - thetas2(addPiTh);
    [thetas2, indS] = sort(thetas2);
    sigma2 = sigma2(indS);
    % Plot
    plot(thetas2(:), sigma2(:), '-o', 'Color', colors(i, :))
end


leg = string(xcoordCtr2 - 1.5);
xlabel("Angle")
ylabel("\sigma")
legendS = legend(leg);
title(legendS, 'd')
title("Solution density and angle, \Omega_{2}")





