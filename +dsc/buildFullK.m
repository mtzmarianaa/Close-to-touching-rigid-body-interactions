function [K_full,  listKs_inv]= buildFullK(ds, kernel)
% Build the full matrix for a BIE using kernel as a kernel.
% K_full = [K11 K12; K21 K22]


Ntot = ds.chnkrs.npt;
nB = ds.nB;
nChunkers = length(ds.listChnkrs);

K_full = zeros(Ntot);
opts = [];
opts2 = [];
opts2.adaptive_correction = true;

if( nargout < 2 )
    % Just build the full K matrix, don't store or compute inverses
     % Proceed to fill in the matrix 
    for k=1:nChunkers
        % fill column by column
        chnkrk = ds.listChnkrs(k);
        targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
        start_col = nB(k) + 1;
        end_col = nB(k+1);
        for i=1:nChunkers
            start_row = nB(i) + 1;
            end_row = nB(i + 1);
            chnkri = ds.listChnkrs(i); % chunker we are working with
            % See if we have to do an off boundary or on boundary eval
            if(i == k)
                % on boundary
                submat = chunkermat(chnkri, kernel, opts2);
            else
                % off boundary
                submat = chunkerkernevalmat(chnkri, kernel, targ, opts);
            end
            % Add this submatrix to K
            K_full(start_col:end_col, start_row:end_row) = submat;
        end
    end   
end

if( nargout == 2 )
    % Compute the full matrix AND store the inverses
    listKs_inv = cell(1, ds.nGammas); % Initialize cell array
     % Proceed to fill in the matrix 
    for k=1:nChunkers
        % fill column by column
        chnkrk = ds.listChnkrs(k);
        targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
        start_col = nB(k) + 1;
        end_col = nB(k+1);
        for i=1:nChunkers
            start_row = nB(i) + 1;
            end_row = nB(i + 1);
            chnkri = ds.listChnkrs(i); % chunker we are working with
            % See if we have to do an off boundary or on boundary eval
            if(i == k)
                % on boundary
                submat = chunkermat(chnkri, kernel);
                if( i > 1)
                    listKs_inv{i-1} = inv(submat); % Save this inverse (it's on the diagonal)
                end
            else
                % off boundary
                submat = chunkerkernevalmat(chnkri, kernel, targ, opts);
            end
            % Add this submatrix to K
            K_full(start_col:end_col, start_row:end_row) = submat;
        end
    end 
end


end