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

[m, n] = size(K_full);

fprintf("\n\n\nNEW GMRES SOLVE FULL \n\n     %5.2e  time solve  \n     " + ...
    "%d x %d matrix dimensions\n     %5.2e condition number\n\n\n", t2, m, n, cond(K_full));

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