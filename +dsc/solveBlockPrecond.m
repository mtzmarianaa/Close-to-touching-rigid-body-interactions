function [sigmaBlockPrecond, nGMRES, tSolve] = solveBlockPrecond( ds, rhs, kernel, matOffSet )
% Solve the BIE for K*sigma = rhs FULLY, with preconditioning but no compressing
% Also outputs the number of GMRES iterations needed to solve the problem

if( nargin < 3)
    matOffSet = sparse(ds.chnkrs.npt, ds.chnkrs.npt);
end
s = tic();
% Build the matrix, inverses are required in this case
[K_full, listKs_inv] = dsc.buildFullK(ds, kernel, matOffSet);
% Build sparse diagonal block with inverses
I0 = eye(ds.gamma0.npt);
block_inv = blkdiag(I0, listKs_inv{:});
K_precond = K_full*block_inv;


[sigmaBlockPrecond, ~, ~, nGMRES] = gmres(K_precond, rhs, [], 1e-10, size(K_precond,1));
t2 = toc(s);
nGMRES = nGMRES(2);
sigmaBlockPrecond = block_inv*sigmaBlockPrecond;

[m,n] = size(K_precond);
fprintf("\n\n\nNEW GMRES SOLVE BLOCK PRECOND \n\n     %5.2e  time solve  \n     " + ...
    "%d x %d matrix dimensions\n     %5.2e condition number\n\n\n", t2, m, n, cond(K_precond));

if(nargout > 2)
    tSolve = t2;
end

end