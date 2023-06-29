function ds = splitDS(ds, id, ich, x, w, u)
% Similar to split in chunkie/chunker. This routine takes the list of all
% chunks in a disc (id) and splits one (ich) in half with respect to the
% parameter space. The new nodes that are added are exactly ON the true
% geometry.
%
% Syntax: ds = split(ds, id, ich, opts, x, w, u)
%
% Input:
%          ds   -   the discs object
%          id    -   the index of the disc where the chunk is
%          ich  -   the chunk number to split
%
% Optional input:
%          x - precomputed Legendre nodes of order chnkr.k
%          w - precomputed Legendre weights
%          u - precomputed vals at legendre nodes -> coeffs matrix
%
% Output:
%          ds  - modified discs object
%
% TODO: maybe consider splitting the chunks in another way? CURRENTLY IN
% PAUSE, USING SPLING FROM CHUNKIE

if( isempty(ds.thetas(id)) )
    fprintf("\nCan't split disc, no information about thetas saved. \n");
    return
end

chnkr = ds.listChnkrs(id); % Get the chunker object we are dealing with

if( nargin < 6)
    % Need to compute Legendre nodes, weights, matrix 
    [x, w, u] = lege.exps(chnkr.k);
end

% Get the parameters we are interested in changing
r = chnkr.r(:, :, ich);
d = chnkr.d(:, :, ich);
d2 = chnkr.d2(:, :, ich);
nch = chnkr.nch;

% 




end