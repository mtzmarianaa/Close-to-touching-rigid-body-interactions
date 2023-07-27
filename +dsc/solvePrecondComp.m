function [sigmaPrecondComp, nGMRES] = solvePrecondComp( ds, rhsC, kernel, matOffSet, matOffSetCoarse, listKs_inv, listP, listR )
% Solve the BIE for K*sigma = rhs FULLY, with preconditioning AND compressing
% Also outputs the number of GMRES iterations needed to solve the problem

% Options for on boundary evaluations
opts2 = [];
opts2.adaptive_correction = true;

nGammas = ds.nGammas;
nChunkers = length(ds.listChnkrs); % Total number of chunks
opts = []; % Options for building off diagonal matrices with chunkie

% Build nBc, nB for the compressed blocks
nBc = ds.nBCoarse;
nB = ds.nB;

if( nargin < 5)
    matOffSet = sparse(ds.chnkrs.npt, ds.chnkrs.npt);
    matOffSetCoarse = sparse( ds.nBCoarse(end), ds.nBCoarse(end) );
end

% Build the inverses and the R matrices if needed
if( nargin < 6 )
    listKs_inv = cell(1, ds.nGammas);
    listP = cell(1, ds.nGammas);
    listR = cell(1, ds.nGammas);
    for i =1:nGammas
        K22 = chunkermat(ds.listGammas(i), kernel, opts2);
        K22_inv = inv(  K22 + matOffSet( (nB(i+1)+1):nB(i+2),  (nB(i+1)+1):nB(i+2) ) );
        listKs_inv{i} = K22_inv;
        nRef = floor(ds.listGammas(i).nch/4 - 2);
        nRef = max(0, nRef);
        P = rcip.prol_dyadic(ds.listCoarseGammas(i).k, nRef);
        P = blkdiag(P, P);
        listP{i} = P;
        listR{i} = rcip.buildR(ds.listCoarseGammas(i), ds.listGammas(i), K22_inv, P, kernel);
    end
end

if( nargin < 8 )
    listR = cell(1, ds.nGammas);
        for i =1:nGammas
            listR{i} = rcip.buildR(ds.listCoarseGammas(i), ds.listGammas(i), listKs_inv{i}, listP{i}, kernel);
        end
end

NtotC = ds.gamma0.npt + sum( ds.listCoarseGammas(:).npt );

% Build the "big" identity matrix
bigEye = speye(NtotC);

% Build the block diagonal matrix
I0 = eye(ds.gamma0.npt);
block_R = blkdiag(I0, listR{:});
blockProlongation = blkdiag(I0, listP{:});

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
[sigmaPrecondComp, ~, ~, nGMRES] = gmres( Kc, rhsC, [], 1e-10, size(Kc,1) );
t2 = toc(s);

[m,n] = size(Kc);
fprintf("\n\n\nNEW GMRES SOLVE BLOCK PRECOND COMPR \n\n     %5.2e  time solve  \n     " + ...
    "%d x %d matrix dimensions\n     %5.2e condition number\n\n\n", t2, m, n, cond(Kc));

nGMRES = nGMRES(2);

% Get sigma in the fine discretization


sigmaPrecondComp = blockProlongation*block_R*sigmaPrecondComp;


end