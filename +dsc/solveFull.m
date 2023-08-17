function [sigmaFull, nGMRES, tSolve] = solveFull( ds, rhs, kernel, matOffSet, verbose)
% *solveFull* solve the BIE for K*sigma = rhs. Solves the full
% problem as stated. Solves the linear system with
% GMRES. Also outputs number of GMRES iterations needed to solve the
% problem. Could output time taken to assemble and solve the problem.
%
% Syntax: [sigmaFull, nGMRES] = solveFull( ds, rhs, kernel, matOffSet )
%              [sigmaFull, nGMRES, tSolve] = solveFull( ds, rhs, kernel, matOffSet )
%              [sigmaFull, nGMRES, tSolve] = solveFull( ds, rhs, kernel, matOffSet, verbose)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   rhs - right hand side of the BIE
%   kernel - kernel object (from chunkie) or function handle defining the
%                kernel to use
%
% Optional input:
%   matOffSet - matrix defining the integral operator which is not a kernel
%                        (has to be a matrix)
%   verbose - true or false, if print reassuring statements
%
% Output:
%   sigmaFull - solution density for the BIE (organized by blocks)
%   nGMRES - number of GMRES iterations needed to solve the linear system
%
% Optional output:
%   tSolve - time taken to assemble and solve the problem
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

if( isa(kernel, 'kernel') )
    kEval = @(s,t) kernel.eval(s,t);
else
    kEval = @(s,t) kernel(s,t);
end


if( nargin < 4)
    matOffSet = zeros(ds.chnkrs.npt);
end

if(nargin<5)
    verbose = true;
end

s = tic();
% Build the matrix, no inverses required
K_full = dsc.buildFullK(ds, kEval, matOffSet);


[sigmaFull, ~, ~, nGMRES] = gmres(K_full, rhs, [], 1e-10, size(K_full, 1));
t2 = toc(s);

nGMRES = nGMRES(2);

[m, n] = size(K_full);
if(verbose)
fprintf("\n\n\nNEW GMRES SOLVE FULL \n\n     %5.2e  time solve  \n     " + ...
    "%d x %d matrix dimensions\n     %5.2e condition number\n\n\n", t2, m, n, cond(K_full));
end

if(nargout > 2)
    tSolve = t2;
end

% cmap_bpp = rgbmap('white', 'powder blue', 'cornflower blue', 'purple blue', 'royal purple', 1024);
% figure()
% x = [0, m];
% y = [0, m];
% imagesc(x, y, K_full)
% title("K_{full}")
% colorbar()
% colormap(cmap_bpp)
% axis image

end