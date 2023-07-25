function [sigmaInterpPrecond, nGMRES] = solveInterpPrecond(ds, rhsC, kernel, matInterpolant, matOffSet, matOffSetCoarse, geom0)
% Solve the BIE for K*sigma = rhs FULLY, with interpolation, preconditioning AND compressing
% Also outputs the number of GMRES iterations needed to solve the problem
% TWO DISCS ONLY

assert( ds.nGammas==1, 'Interpolation, preconditioning and compression only supported for two discs at this moment \n' );

if( nargin < 5 )
    matOffSetCoarse = zeros( (ds.gamma0.npt + sum(ds.listCoarseGammas(:).npt)) );
end

if( nargin < 7 )
    geom0 = [];
    geom0.Rs = [0.75; 0.75];
    geom0.ctrs = [0  1.6; 0 0];
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

if( nargin < 6)
    matOffSet = sparse(ds.chnkrs.npt, ds.chnkrs.npt);
    matOffSetCoarse = sparse( ds.nBCoarse(end), ds.nBCoarse(end) );
end


% Interpolation part
d = norm( ds.ctrs(:, 2) - ds.ctrs(:, 1)) - ds.Rs(1) - ds.Rs(2);
kCh = floor( log( (norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2))/d  )*(1/log(2)) );
%   corresponding coefficients

da = 1/(2^kCh)*(norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2) );
db = 1/(2^(kCh + 1))*(norm(geom0.ctrs(:, 1) - geom0.ctrs(:, 2)) - geom0.Rs(1) - geom0.Rs(2) );

xDist = 2*(d - da)/(db - da) - 1;

% Evaluate the interpolated R
R_interpolated = rcip.evaluateRInterpolated(xDist, matInterpolant{kCh+1});


% Build the system
nRef = floor(ds.listGammas(1).nch/4 - 2);
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
K11 = chunkermat( ds.gamma0, kernel, opts2) + matOffSet(1:nB(2), 1:nB(2));
K11 = K11 - eye(size(K11));
Kc(1:nB(2), 1:nB(2)) = K11;
% Build the first row and column of the Kc matrix
for i=2:nChunkers
    targ = reshape( ds.listCoarseGammas(i-1).r, 2, []);
    rowBlock = chunkerkernevalmat( ds.gamma0, kernel, targ, opts  );
    Kc(  (nBc(i)+1):nBc(i+1), 1:nB(2)  ) = rowBlock + matOffSetCoarse(  (nBc(i)+1):nBc(i+1), 1:nB(2)  );
    targGamma0 = reshape( ds.gamma0.r, 2, [] );
    columnBlock = chunkerkernevalmat( ds.listCoarseGammas(i-1), kernel, targGamma0, opts );
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
            Keval = chunkerkernevalmat(chnkri, kernel, targ, opts);
            Kc(start_col:end_col, start_row:end_row) = Keval + matOffSetCoarse(start_col:end_col, start_row:end_row);
        end
    end
end   
%%%%%%%%%%%%%%%%%%%%%%%%

% Put the 3 matrices together
Kc = bigEye + Kc*block_R;

% Solve the linear system using GMRES

s = tic();
[sigmaInterpPrecond, ~, ~, nGMRES] = gmres( Kc, rhsC, [], 1e-10, size(Kc,1) );
t2 = toc(s);
fprintf("%5.2e s :time taken to solve the compressed precond linear system with GMRES\n", t2);
nGMRES = nGMRES(2);

% Get sigma in the fine discretization
blockProlongation = blkdiag(I0, P);
sigmaInterpPrecond = blockProlongation*block_R*sigmaInterpPrecond;


end