% Generates plots comparing different methods: full, preconditioning,
% compression + preconditioning, interpolation + preconditioning +
% compression

nTest = 8;

tSolve_capacitanceFull = zeros(nTest, 1);
tSolve_elastanceFull = zeros(nTest, 1);
tSolve_capacitancePrecond = zeros(nTest, 1);
tSolve_elastancePrecond = zeros(nTest, 1);
tSolve_capacitancePrecondCompress = zeros(nTest, 1);
tSolve_elastancePrecondCompress = zeros(nTest, 1);
tSolve_capacitanceInterPrecondCompress = zeros(nTest, 1);
tSolve_elastanceInterPrecondCompress = zeros(nTest, 1);

% Full
for i=1:10
    CapacitanceElastanceErrorsAdaptive
    d = xcoordCtr2-1.5;
    nGMRESFull = nGMRES_capacitance;
    errorsFull = errors_ukAdaptive;
    tSolve_capacitanceFull = tSolve_capacitanceFull + tSolve_capacitance;
    tSolve_elastanceFull = tSolve_elastanceFull + tSolve_elastance;
    close all
end
tSolve_capacitanceFull = tSolve_capacitanceFull./10;
tSolve_elastanceFull = tSolve_elastanceFull./10;

% Preconditioning
for i=1:10
    CapacitanceElastanceErrorsAdaptivePrecond
    nGMRESPrecond = nGMRES_capacitance;
    errorsPrecond = errors_ukAdaptive;
    tSolve_capacitancePrecond = tSolve_capacitancePrecond + tSolve_capacitance;
    tSolve_elastancePrecond = tSolve_elastancePrecond + tSolve_elastance;
    close all
end
tSolve_capacitancePrecond = tSolve_capacitancePrecond./10;
tSolve_elastancePrecond = tSolve_elastancePrecond./10;

% Compression + Preconditioning
for i=1:10
    CapacitanceElastanceErrorsAdaptiveCompressPrecond
    nGMRESPrecondCompress = nGMRES_capacitance;
    errorsPrecondCompress = errors_ukAdaptive;
    tSolve_capacitancePrecondCompress = tSolve_capacitancePrecondCompress + tSolve_capacitance;
    tSolve_elastancePrecondCompress = tSolve_elastancePrecondCompress + tSolve_elastance;
    close all
end
tSolve_capacitancePrecondCompress = tSolve_capacitancePrecondCompress./10;
tSolve_elastancePrecondCompress = tSolve_elastancePrecondCompress./10;

% Interpolation + Compression + Preconditioning
for i=1:10
    CapacitanceElastanceErrorsAdaptiveInterpCompressPrecond
    nGMRESInterPrecondCompress = nGMRES_capacitance;
    errorsInterPrecondCompress = errors_ukAdaptive;
    tSolve_capacitanceInterPrecondCompress = tSolve_capacitanceInterPrecondCompress + tSolve_capacitance;
    tSolve_elastanceInterPrecondCompress = tSolve_elastanceInterPrecondCompress + tSolve_elastance;
    close all
end
tSolve_capacitanceInterPrecondCompress = tSolve_capacitanceInterPrecondCompress./10;
tSolve_elastanceInterPrecondCompress = tSolve_elastanceInterPrecondCompress./10;

% Plot comparisons

colors = [0 0 160; 0 155 255; 120 0 255; 0 255 255];
colors = colors./255;

% nGMRES
fig1 = figure()
loglog(d, nGMRESFull, '-o', 'Color', colors(1, :))
hold on
loglog(d, nGMRESPrecond, '-o', 'Color', colors(2, :))
loglog(d, nGMRESPrecondCompress, '-o', 'Color', colors(3, :))
loglog(d, nGMRESInterPrecondCompress, '-o', 'Color', colors(4, :))
title('Distance between discs, GMRES iterations')
xlabel('d')
ylabel('nGMRES')
hold off
leg = ["Full", "Block Preconditioning", "Preconditioning and Compression", "Interpolation, Preconditioning, and Compression"];
legendS = legend(leg);
title(legendS, 'Solver used')
xlim([0, max(d)])
fontsize(fig1, 24, "points")


% Errors
fig2 = figure()
loglog(d, errorsFull, '-o', 'Color', colors(1, :))
hold on
loglog(d, errorsPrecond, '-o', 'Color', colors(2, :))
loglog(d, errorsPrecondCompress, '-o', 'Color', colors(3, :))
loglog(d, errorsInterPrecondCompress, '-o', 'Color', colors(4, :))
title('Distance between discs, errors_{rel}')
xlabel('d')
ylabel('errors_{rel}')
leg = ["Full", "Block Preconditioning", "Preconditioning and Compression", "Interpolation, Preconditioning, and Compression"];
legendS = legend(leg);
title(legendS, 'Solver used')
xlim([0, max(d)])
fontsize(fig2, 24, "points")

% tSolve
fig3 = figure()
loglog(d, tSolve_capacitanceFull, '-o', 'Color', colors(1, :))
hold on
loglog(d, tSolve_capacitancePrecond, '-o', 'Color', colors(2, :))
loglog(d, tSolve_capacitancePrecondCompress, '-o', 'Color', colors(3, :))
loglog(d, tSolve_capacitanceInterPrecondCompress, '-o', 'Color', colors(4, :))
title('Distance between discs, time')
xlabel('d')
ylabel('T')
leg = ["Full", "Block Preconditioning", "Preconditioning and Compression", "Interpolation, Preconditioning, and Compression"];
legendS = legend(leg);
title(legendS, 'Solver used')
xlim([0, max(d)])
fontsize(fig3, 24, "points")
ylim([0, max(tSolve_capacitanceFull) + 5])




