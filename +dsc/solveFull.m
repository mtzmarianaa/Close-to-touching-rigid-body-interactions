function [sigmaFull, nGMRES] = solveFull( ds, rhs, kernel )
% Solve the BIE for K*sigma = rhs FULLY, no preconditioning or compressing
% Also outputs the number of GMRES iterations needed to solve the problem

% Build the matrix, no inverses required
K_full = dsc.buildFullK(ds, kernel);

s = tic();
[sigmaFull, nGMRES] = gmres(K_full, rhs, [], 1e-14, 1500);
t2 = toc(s);

fprintf("%5.2e s :time taken to solve the full linear system with GMRES\n", t2);

end