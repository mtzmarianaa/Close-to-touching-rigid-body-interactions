function [K_full,  listKs_inv]= buildFullK(ds, kernel,  matOffSet)
% *buildFullK* builds the full matrix for a BIE using kernel as a kernel.
% Has the option to set up a matOffSet (integral operator not defined by a
% "standard" kernel. Can also output the list of inverses of the
% close-to-touching regions (useful for RCIP).
%
% Syntax: [K_full,  listKs_inv]= buildFullK(ds, kernel)
%              [K_full,  listKs_inv]= buildFullK(ds, kernel,  matOffSet)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   kernel - kernel object (from chunkie) or function handle definind the
%                kernel to use
%
% Optional input:
%   matOffSet - matrix defining the integral operator which is not a kernel
%                        (has to be a matrix)
%
% Output:
%   K_full - (ds.chnkrs.npt , ds.chnkrs.npt) block matrix (first row and
%                column are the interactions with the far regions of the discs
%                $$\begin{bmatrix} K_{11} & K_{12} & \hdots & K_{1, N} \\ 
%                                               K_{21} & K_{22} & \hdots & K_{2, N} \\ 
%                                               \vdots & \vdots & \ddots & \vdots \\ 
%                                                K_{N1} & K_{N2} & \hdots & K_{N, N} \\ \end{bmatrix} $$
% Optional output:
%   listKs_inv - list of inverses for the close-to-touching region,
%                      $K_{22}^{-1}, K_{33}^{-1}, ..., K_{NN}^{-1}$
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

% Determine if matrix is given, if not initialize it with zeros
if( nargin < 3)
    matOffSet = sparse(ds.chnkrs.npt, ds.chnkrs.npt);
end


Ntot = ds.chnkrs.npt;
nB = ds.nB;
nChunkers = length(ds.listChnkrs);

K_full = zeros(Ntot);
opts = [];
opts2 = [];
opts2.adaptive_correction = true;

% If we don't need the inverses just fill the matrix block by block
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

% If requested, compute the inverses for the close-to-touching regions
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
                submat = chunkermat(chnkri, kernel, opts2);
                if( i > 1)
                    listKs_inv{i-1} = inv(submat + matOffSet(start_col:end_col, start_row:end_row)); % Save this inverse (it's on the diagonal)
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

% Assemble the whole matrix
K_full = K_full + matOffSet;

end