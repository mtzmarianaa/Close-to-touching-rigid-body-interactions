function [chnkr1, chnkr2] = divideChnkr(chnkr, opts, pref)

% Divides a chnkr object into two smaller chnkrs. 
%
% Inputs:
% chnkr - the chnkr object to divide
%
% Optional inputs:
% indices - int(2), the indices on which to divide the chunkers (ceil(chnkr.nch)/2)



if(nargin < 2)
    opts = [];
    opts.indices = ceil(chnkr.nch/2);
end


if(nargin < 3)
    pref = chunkerpref();
end

% Get the number of chunks in each chunker
nch = chnkr.nch;
nch1 = opts.indices;
nch2 = nch - nch1;


% Initialize the two chunkers and add their corresponding number of chunks
chnkr1 = chunker(pref);
chnkr1 = chnkr1.addchunk(nch1);
chnkr2 = chunker(pref);
chnkr2 = chnkr2.addchunk(nch2);

% Now set the information accordingly BE CAREFUL WITH ADJ

chnkr1.r = chnkr.r(:, :, 1:nch1);
chnkr1.d = chnkr.d(:, :, 1:nch1);
chnkr1.d2 = chnkr.d2(:, :, 1:nch1);
chnkr1.h = chnkr.h(1:nch1);
chnkr1.adj = chnkr.adj(:, 1:nch1);
% Change the adjacent panels
chnkr1.adj(1, 1) = chnkr1.adj(2, 1);
chnkr1.adj(2, nch1) = chnkr1.adj(1, nch1);

chnkr2.r = chnkr.r(:, :, (nch1+1):nch);
chnkr2.d = chnkr.d(:, :, (nch1+1):nch);
chnkr2.d2 = chnkr.d2(:, :, (nch1+1):nch);
chnkr2.h = chnkr.h( (nch1+1):nch );
chnkr2.adj = chnkr.adj(:, (nch1+1):nch );
% Change the adjacent panels
chnkr2.adj(2, 1) = chnkr2.adj(1, 1);
chnkr2.adj = chnkr2.adj - nch2;
chnkr2.adj(1, nch2) = nch2 - 1;
chnkr2.adj(2, nch2) = nch2 - 1;



end