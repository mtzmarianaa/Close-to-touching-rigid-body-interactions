function [sigmaFull, nGMRES] = solveFull( ds, rhs, kernel, matOffSet )
% Solve the BIE for K*sigma = rhs FULLY, no preconditioning or compressing
% Also outputs the number of GMRES iterations needed to solve the problem


if( nargin < 3)
    matOffSet = zeros(ds.chnkrs.npt);
end

% Build the matrix, no inverses required
K_full = dsc.buildFullK(ds, kernel, matOffSet);

s = tic();
[sigmaFull, ~, ~, nGMRES] = gmres(K_full, rhs, [], 1e-10, size(K_full, 1));
t2 = toc(s);

nGMRES = nGMRES(2);

fprintf("%5.2e s :time taken to solve the full linear system with GMRES\n", t2);

end