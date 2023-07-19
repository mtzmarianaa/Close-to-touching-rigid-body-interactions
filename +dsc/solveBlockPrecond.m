function [sigmaBlockPrecond, nGMRES] = solveBlockPrecond( ds, rhs, kernel, matOffSet )
% Solve the BIE for K*sigma = rhs FULLY, with preconditioning but no compressing
% Also outputs the number of GMRES iterations needed to solve the problem

if( nargin < 3)
    matOffSet = sparse(ds.chnkrs.npt, ds.chnkrs.npt);
end

% Build the matrix, inverses are required in this case
[K_full, listKs_inv] = dsc.buildFullK(ds, kernel, matOffSet);
% Build sparse diagonal block with inverses
I0 = eye(ds.gamma0.npt);
block_inv = blkdiag(I0, listKs_inv{:});
K_precond = K_full*block_inv;

s = tic();
[sigmaBlockPrecond, ~, ~, nGMRES] = gmres(K_precond, rhs, [], 1e-10, size(K_precond,1));
t2 = toc(s);
nGMRES = nGMRES(2);
sigmaBlockPrecond = block_inv*sigmaBlockPrecond;

fprintf("%5.2e s :time taken to solve the full precond linear system with GMRES\n", t2);

end