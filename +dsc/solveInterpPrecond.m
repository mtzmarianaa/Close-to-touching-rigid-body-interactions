function [sigmaInterpPrecond, nGMRES, tSolve] = solveInterpPrecond(ds, rhsC, kernel, listPrecomputedR, typeNodes, matOffSet, matOffSetCoarse, geom0, verbose)
% *solveInterpPrecond* solve the BIE for K*sigma = rhs. Solves the compressed,
% preconditioned system using interpolation from precomputed matrices. 
% Solves the linear system with GMRES. Also outputs number of GMRES 
% iterations needed to solve the problem. 
% Could output time taken to assemble and solve the problem.
%
% Syntax: [sigmaInterpPrecond, nGMRES, tSolve] = solveInterpPrecond(ds, rhsC, kernel, listPrecomputedR, typeNodes)
%              [sigmaInterpPrecond, nGMRES, tSolve] = solveInterpPrecond(ds, rhsC, kernel, listPrecomputedR, typeNodes, matOffSet, matOffSetCoarse)
%              [sigmaInterpPrecond, nGMRES, tSolve] = solveInterpPrecond(ds, rhsC, kernel, listPrecomputedR, typeNodes, matOffSet, matOffSetCoarse, geom0)
%              [sigmaInterpPrecond, nGMRES, tSolve] = solveInterpPrecond(ds, rhsC, kernel, listPrecomputedR, typeNodes, matOffSet, matOffSetCoarse, geom0, verbose)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   rhs - right hand side of the BIE
%   kernel - kernel object (from chunkie) or function handle definind the
%                kernel to use
%   listPrecomputedR - list of precomputed R matrices (has to be given)
%   typeNodes - 'l' for Legendre nodes, 'logc' for log Chebyshev
%
% Optional input:
%   matOffSet - matrix defining the integral operator in the fine mesh which is not a kernel
%                        (has to be a matrix)
%   matOffSetCoarse - matrix defining the integral operator in the coarse mesh which is not a kernel
%                        (has to be a matrix)
%   geom0 - initial geometry of the discs for distance0 (for the
%                interpolation on the distance between two discs)
%   verbose - true or false, if print reassuring statements
%
% Output:
%   sigmaInterpPrecond - solution density for the BIE (organized by blocks)
%   nGMRES - number of GMRES iterations needed to solve the linear system
%
% Optional output:
%   tSolve - time taken to assemble and solve the problem
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

%assert( ds.nGammas==1, 'Interpolation, preconditioning and compression only supported for two discs at this moment \n' );

if( isa(kernel, 'kernel') )
    kEval = @(s,t) kernel.eval(s,t);
else
    kEval = @(s,t) kernel(s,t);
end


% Determine if matrix is given, if not initialize it with zeros
if( nargin < 7 || length(matOffSet(:)) ==1 || length(matOffSetCoarse(:)) == 1) 
    matOffSet = sparse(ds.chnkrs.npt, ds.chnkrs.npt);
    matOffSetCoarse = sparse( ds.nBCoarse(end), ds.nBCoarse(end) );
end

% Determine if initial geometry is given, otherwise just use what we have
if( nargin < 8 || ~isstruct(geom0) )
    geom0 = [];
    geom0.Rs = [0.75; 0.75];
    geom0.ctrs = [0  1.6; 0 0];
end

if(nargin < 9)
    verbose = true;
end

% Options for on boundary evaluations
opts2 = [];
opts2.adaptive_correction = true;

nGammas = ds.nGammas;
nChunkers = length(ds.listChnkrs); % Total number of chunks
opts = []; % Options for building off diagonal matrices with chunkie

% Build nBc, nB for the compressed blocks
nBc = ds.nBCoarse;
nB = ds.nB;

s = tic();
% Interpolation part
d = norm( ds.ctrs(:, 2) - ds.ctrs(:, 1)) - ds.Rs(1) - ds.Rs(2);

if strcmp(typeNodes, 'l')
    kCh = floor( log( (norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2))/d  )*(1/log(2)) );
    kCh = max(0, kCh);
    %   corresponding coefficients
    da = 1/(2^kCh)*(norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2) );
    db = 1/(2^(kCh + 1))*(norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2) );
    xDist = 2*(d - da)/(db - da) - 1;
else
    xDist = -1 + 2*(d - 0.1)/(1e-12 - 0.1);
    kCh = 0;
end


% Evaluate the interpolated R
R_interpolated = rcip.evaluateRInterpolated(xDist, listPrecomputedR{kCh+1}, typeNodes);

% Build the system
nRef = floor(ds.listGammas(1).nch/4 - 2);
nRef = max(0, nRef);
P = rcip.prol_dyadic(ds.listCoarseGammas(1).k, nRef);
P = blkdiag(P, P);
NtotC = ds.gamma0.npt + sum( ds.listCoarseGammas(:).npt );

% Build the "big" identity matrix
bigEye = speye(NtotC);

% Build the block diagonal matrix
I0 = eye(ds.gamma0.npt);
block_R = blkdiag(I0, R_interpolated);

%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kc matrix
Kc = zeros( NtotC );

% First block
K11 = chunkermat( ds.gamma0, kEval, opts2) + matOffSet(1:nB(2), 1:nB(2));
K11 = K11 - eye(size(K11));
Kc(1:nB(2), 1:nB(2)) = K11;
% Build the first row and column of the Kc matrix
for i=2:nChunkers
    targ = reshape( ds.listCoarseGammas(i-1).r, 2, []);
    rowBlock = chunkerkernevalmat( ds.gamma0, kEval, targ, opts  );
    Kc(  (nBc(i)+1):nBc(i+1), 1:nB(2)  ) = rowBlock + matOffSetCoarse(  (nBc(i)+1):nBc(i+1), 1:nB(2)  );
    targGamma0 = reshape( ds.gamma0.r, 2, [] );
    columnBlock = chunkerkernevalmat( ds.listCoarseGammas(i-1), kEval, targGamma0, opts );
    Kc( 1:nB(2), (nBc(i)+1):nBc(i+1) ) = columnBlock + matOffSetCoarse( 1:nB(2), (nBc(i)+1):nBc(i+1) );
end

% Build the second block - the one with blocks of zeros on the diagonal
for k=1:nGammas
    % fill column by column
    chnkrk = ds.listCoarseGammas(k);
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nB(k+1) + 1;
    end_col = nB(k+2);
    for i=1:nGammas
        start_row = nB(i+1) + 1;
        end_row = nB(i + 2);
        chnkri = ds.listCoarseGammas(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i ~= k)
            % off boundary, else its just zeros and we don't need to
            % compute anything
            Keval = chunkerkernevalmat(chnkri, kEval, targ, opts);
            Kc(start_col:end_col, start_row:end_row) = Keval + matOffSetCoarse(start_col:end_col, start_row:end_row);
        end
    end
end

% Put the 3 matrices together
Kc = bigEye + Kc*block_R;

% Solve the linear system using GMRES
[sigmaInterpPrecond, ~, ~, nGMRES] = gmres( Kc, rhsC, [], 1e-10, size(Kc,1) );
t2 = toc(s);

% cmap_bpp = rgbmap('white', 'powder blue', 'cornflower blue', 'purple blue', 'royal purple');
% figure()
% imagesc(inv(block_R)*Kc)
% title("K_interp")
% colorbar()
% colormap(cmap_bpp)

[m,n] = size(Kc);
if(verbose)
fprintf("\n\n\nNEW GMRES SOLVE BLOCK PRECOND INTERPOLATION \n\n     %5.2e  time solve  \n     " + ...
    "%d x %d matrix dimensions\n     %5.2e condition number\n\n\n", t2, m, n, cond(Kc));
end
nGMRES = nGMRES(2);

% Get sigma in the fine discretization
blockProlongation = blkdiag(I0, P);
sigmaInterpPrecond = blockProlongation*block_R*sigmaInterpPrecond;

if(nargout > 2)
    tSolve = t2;
end

end