%%% Compare the errors as discs come together. Solve first capacitance,
%%% then elastance and compare
% FULL SYSTEM SOLVE
addpaths_loc();
% clear all
% close all
clc

u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);

solveType = 'full';

uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 1.5 ;0 0]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;


%%%%%%%%%%%%%%%%%%%%
%%%% Adaptive

pClose = [];
pClose(1).data = [0 2 1];
pClose(1).nClose = 1;
pClose(1).thetasReg = pi/6;
pClose(2).data = [pi, 1, 1];
pClose(2).nClose =1;
pClose(2).thetasReg = pi/6;


xcoordCtr2 = [ 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8  ];
nTest = length(xcoordCtr2);


errors_ukAdaptive = zeros( nTest, 1 );
nGMRES_capacitance = zeros( nTest, 1);
nGMRES_elastance = zeros( nTest, 1 );
depth = zeros(nTest, 1);

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    geom.ctrs = ctrs;
    ds = discs(geom, pClose);
    depth(i) = floor( ds.listGammas(1).nch/4 - 2  );
    % First solve capacitance
    [qkC, sigmaC, nGMRES_C] = capacitanceProblem(ds, uk, solveType);
    % Compute ukC
    ukC = zeros(ds.chnkrs.npt, 1);
    for k = 1:n 
        % Use the flag
        flag = logical( dsc.flagnDisc(k, ds) );
        % Convert that to points
        flag_points = repmat(flag, 1, ds.chnkrs.k);
        flag_points = flag_points';
        flag_points = logical( flag_points(:) );
        % Get the points
        xOnSurface = reshape(ds.chnkrs.r(:, :, flag ), 2, []);
        u_toUse = uk{k}(xOnSurface);
        ukC(flag_points) = u_toUse';
    end
    % Now with this q solve elastance
    [ukE, sigmaE, nGMRES_E] = elastanceProblem(ds, qkC, solveType);
    % Compute the difference
    errors_ukAdaptive(i) = norm(ukC - ukE)/norm(ukC);
    % Add the number of iterations needed
    nGMRES_capacitance(i) = nGMRES_C;
    nGMRES_elastance(i) = nGMRES_E;
end



% Plot
cq = [0 232/255 255/255];
cO = [147/255 155/255 255/255];
cS = [155/255 0 255/255];

figure()
plot(xcoordCtr2-1.5, errors_ukAdaptive, '-o', 'Color', cq)
title("Errors on surface density - Full")
xlabel("Distance between discs")
ylabel("Relative error (l2)")

figure()
loglog(xcoordCtr2-1.5, errors_ukAdaptive, '-o', 'Color', cq)
title("Errors on surface density (log log) - Full")
xlabel("Distance between discs")
ylabel("Relative error (l2)")

figure()
loglog(depth, errors_ukAdaptive, '-o', 'Color', cq)
title("Errors on surface density (log log) - Full")
xlabel("Depth of refinement")
ylabel("Relative error (l2)")


figure()
plot(xcoordCtr2-1.5, nGMRES_capacitance, '-*', 'Color', cO)
title("Number of GMRES iterations needed - Capacitance - Full")
xlabel("Distance between discs")
ylabel("GMRES iterations needed")

figure()
plot(depth, nGMRES_capacitance, '-*', 'Color', cO)
title("Number of GMRES iterations needed - Capacitance - Full")
xlabel("Depth of refinement")
ylabel("GMRES iterations needed")


figure()
plot(xcoordCtr2-1.5, nGMRES_elastance, '-*', 'Color', cS)
title("Number of GMRES iterations needed - Elastance - Full")
xlabel("Distance between discs")
ylabel("GMRES iterations needed")


figure()
plot(depth, nGMRES_elastance, '-*', 'Color', cS)
title("Number of GMRES iterations needed - Elastance - Full")
xlabel("Depth of refinement")
ylabel("GMRES iterations needed")



