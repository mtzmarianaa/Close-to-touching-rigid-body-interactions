% Generates plots comparing different methods: full, preconditioning,
% compression + preconditioning, interpolation + preconditioning +
% compression

% Full
CapacitanceElastanceErrorsAdaptive
d = xcoordCtr2-1.5;
nGMRESFull = nGMRES_capacitance;
errorsFull = errors_ukAdaptive;

% Preconditioning
CapacitanceElastanceErrorsAdaptivePrecond
nGMRESPrecond = nGMRES_capacitance;
errorsPrecond = errors_ukAdaptive;

% Compression + Preconditioning
CapacitanceElastanceErrorsAdaptiveCompressPrecond
nGMRESPrecondCompress = nGMRES_capacitance;
errorsPrecondCompress = errors_ukAdaptive;

% Interpolation + Compression + Preconditioning
CapacitanceElastanceErrorsAdaptiveInterpCompressPrecond
nGMRESInterPrecondCompress = nGMRES_capacitance;
errorsInterPrecondCompress = errors_ukAdaptive;


% Plot comparisons

colors = [0 0 160; 0 155 255; 120 0 255; 0 255 255];
colors = colors./255;

figure()
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


figure()
loglog(d, errorsFull, '-o', 'Color', colors(1, :))
hold on
loglog(d, errorsPrecond, '-o', 'Color', colors(2, :))
loglog(d, errorsPrecondCompress, '-o', 'Color', colors(3, :))
loglog(d, errorsInterPrecondCompress, '-o', 'Color', colors(4, :))
title('Distance between discs, relative errors')
xlabel('d')
ylabel('nGMRES')
leg = ["Full", "Block Preconditioning", "Preconditioning and Compression", "Interpolation, Preconditioning, and Compression"];
legendS = legend(leg);
title(legendS, 'Solver used')
xlim([0, max(d)])