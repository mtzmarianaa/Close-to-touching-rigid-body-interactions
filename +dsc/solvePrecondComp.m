function [sigmaPrecondComp, nGMRES] = solvePrecondComp( ds, rhsC, kernel, listKs_inv, listP, listR )
% Solve the BIE for K*sigma = rhs FULLY, with preconditioning AND compressing
% Also outputs the number of GMRES iterations needed to solve the problem

nGammas = ds.nGammas;
nChunkers = length(ds.listChnkrs); % Total number of chunks
opts = []; % Options for building off diagonal matrices with chunkie

% Build nBc, nB for the compressed blocks
nBc = zeros(nGammas + 2, 1);
nBc(1) = 0;
nBc(2) = ds.gamma0.npt;
for i=3:(nGammas+2)
    nBc(i) = nBc(i-1) + ds.listCoarseGammas(i-2).npt;
end

% Build the inverses and the R matrices if needed
if( nargin < 4 )
    listKs_inv = cell(1, ds.nGammas);
    listP = cell(1, ds.nGammas);
    listR = cell(1, ds.nGammas);
    for i =1:nGammas
        K22_inv = inv(  chunkermat(ds.listGammas(i), kernel) );
        listKs_inv{i} = K22_inv;
        x = lege.exps( ds.listCoarseGammas(i).k );
        nRef = floor(ds.listGammas(i).nch/2);
        Ph = rcip.prol_dyadic(x, nRef);
        P = [Ph Ph];
        listP{i} = P;
        listR{i} = rcip.buildR(ds.listCoarseGammas(i), ds.listGammas(i), kernel, K22_inv, P);
    end
end

if( nargin < 6 )
    listR = cell(1, ds.nGammas);
        for i =1:nGammas
            listR{i} = rcip.buildR(ds.listCoarseGammas(i), ds.listGammas(i), kernel, listKs_inv{i}, listP{i});
        end
end

NtotC = ds.gamma0.npt + sum( ds.listCoarseGammas(:).npt );

% Build the "big" identity matrix
bigEye = speye(NtotC);

% Build the block diagonal matrix
I0 = eye(ds.gamma0.npt);
block_R = blkdiag(I0, listR{:});

%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kc matrix
Kc = zeros( NtotC );

% First block
Kc(1:nBc(1), 1:nBc(1)) = chunkermat( ds.gamma0, kernel ) + eye(nBc(1));
% Build the first row and column of the Kc matrix
for i=2:nChunkers
    targ = reshape( ds.listCoarseGammas(i-1).r, 2, []);
    rowBlock = chunkerkernevalmat( ds.gamma0, kernel, targ, opts  );
    Kc( 1:nBc(1), (nBc(i)+1):nBc(i) ) = rowBlock;
    targGamma0 = reshape( ds.gamma0.r, 2, [] );
    columnBlock = chunkerkernevalmat( ds.listCoarseGammas(i-1), kernel, targGamma0, opts );
    Kc( (nBc(i)+1):nBc(i), 1:nBc(1)  ) = columnBlock;
end

% Build the second block - the one with blocks of zeros on the diagonal
for k=1:nGammas
    % fill column by column
    chnkrk = ds.listCoarseGammas(k);
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nB(k) + 1;
    end_col = nB(k+1);
    for i=1:nGammas
        start_row = nB(i) + 1;
        end_row = nB(i + 1);
        chnkri = ds.listCoarseGammas(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i ~= k)
            % off boundary, else its just zeros and we don't need to
            % compute anything
            Kc(start_col:end_col, start_row:end_row) = chunkerkernevalmat(chnkri, kernel, targ, opts);
        end
    end
end   
%%%%%%%%%%%%%%%%%%%%%%%%

% Put the 3 matrices together

Kc = bigEye + Kc*block_R;

% Solve the linear system using GMRES

s = tic();
[sigmaPrecondComp, nGMRES] = gmres( Kc, rhsC, [], 1e-14 );
t2 = toc(s);
fprintf("%5.2e s :time taken to solve the compressed precond linear system with GMRES\n", t2);


end