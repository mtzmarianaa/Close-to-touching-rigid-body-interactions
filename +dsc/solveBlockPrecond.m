function [sigmaBlockPrecond, nGMRES, tSolve] = solveBlockPrecond( ds, rhs, kernel, matOffSet )
% *solveBlockPrecond* solve the BIE for K*sigma = rhs. Solves the full
% problem using a right preconditioner. Solves the linear system with
% GMRES. Also outputs number of GMRES iterations needed to solve the
% problem. Could output time taken to assemble and solve the problem.
%
% Syntax: [sigmaBlockPrecond, nGMRES] = solveBlockPrecond( ds, rhs, kernel, matOffSet )
%              [sigmaBlockPrecond, nGMRES, tSolve] = solveBlockPrecond( ds, rhs, kernel, matOffSet )
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   rhs - right hand side of the BIE
%   kernel - kernel object (from chunkie) or function handle definind the
%                kernel to use
%
% Optional input:
%   matOffSet - matrix defining the integral operator which is not a kernel
%                        (has to be a matrix)
%
% Output:
%   sigmaBlockPrecond - solution density for the BIE (organized by blocks)
%   nGMRES - number of GMRES iterations needed to solve the linear system
%
% Optional output:
%   tSolve - time taken to assemble and solve the problem
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

% Determine if matrix is given, if not initialize it with zeros
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