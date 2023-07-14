%%% Compare the errors as discs come together. Solve first capacitance,
%%% then elastance and compare
addpaths_loc();
clear all
close all
clc

u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);
nTest = 50;

uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 1.500000000001 ;0 0]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];

%%%%%%%%%%%%%%%%%%%%
%%%% Uniform panels

geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;


%xcoordCtr2 = linspace(1.500000000001, 1.65, nTest ); % Variation in the coordinate ctrs(1,2), moving the discs closer
xcoordCtr2 = 1.5 + [ 1e-3 1e-4 1e-5 1e-6 1e-7  ];
nTest = 5;


errors_uk = zeros( nTest, 1 );
nGMRES_capacitance = zeros( nTest, 1);
nGMRES_elastance = zeros( nTest, 1 );

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    geom.ctrs = ctrs;
    ds = discs(geom);
    % First solve capacitance
    [qkC, sigmaC, nGMRES_C] = capacitanceProblem(ds, uk);
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
    [ukE, sigmaE, nGMRES_E] = elastanceProblem(ds, qkC);
    % Compute the difference
    errors_uk(i) = norm(ukC - ukE)/norm(ukC);
    % Add the number of iterations needed
    nGMRES_capacitance(i) = nGMRES_C(2);
    nGMRES_elastance(i) = nGMRES_E(2);
end



% Plot
cq = [0 232/255 255/255];
cO = [147/255 155/255 255/255];
cS = [155/255 0 255/255];

figure()
plot(xcoordCtr2-1.5, errors_uk, '-o', 'Color', cq)
title("Errors on surface density")
xlabel("Distance between discs")
ylabel("Relative error (l2)")

figure()
loglog(xcoordCtr2-1.5, errors_uk, '-o', 'Color', cq)
title("Errors on surface density (log log)")
xlabel("Distance between discs")
ylabel("Relative error (l2)")


figure()
plot(xcoordCtr2-1.5, nGMRES_capacitance, '-*', 'Color', cO)
title("Number of GMRES iterations needed - Capacitance")
xlabel("Distance between discs")
ylabel("GMRES iterations needed")


figure()
plot(xcoordCtr2-1.5, nGMRES_elastance, '-*', 'Color', cS)
title("Number of GMRES iterations needed - Elastance")
xlabel("Distance between discs")
ylabel("GMRES iterations needed")








