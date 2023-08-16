function [flag] = flagnDiscCoarse(i, ds)
% *flagnDiscCoarse* given a discs object outputs a sparse boolean vector such
% that it indicates which chunk on ds is on the i-th disc. Works only on
% the coarse mesh. See flagnDisc if interested on the fine mesh.
%
% Syntax: [flag] = flagnDiscCoarse(i, ds)
%
% Input:
%   i - index of discs we care about
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%
% Output:
%   flag - (ds.gamma0.nch + sum( ds.listCoarseGammas(:).nch ); , 1) 
%            sparse boolean matrix, 0 if that chunk is
%            not on disc i, 1 if it is on disc i
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

assert(ds.infoClose, 'No information about a coarse and fine discretization of the discs');


nchTot = ds.gamma0.nch + sum( ds.listCoarseGammas(:).nch );

flag = zeros(nchTot, 1); % Initialize the flag

% First part: add those indices in gamma0 that belong to disc i

k = 1;
nBk = zeros( i + 1, 1 );
nBk(1) = 0;
while(k < (i+1) )
    nBk(k+1) = nBk(k) + ds.listFarChunks(k).nch;
    k = k+1;
end
% Fill
flag( nBk(i)+1:nBk(i+1), 1 ) = 1;

% Build indGammasCoarse
indGammasCoarse = zeros(ds.nGammas + 2, 1);
indGammasCoarse(2) = ds.gamma0.nch;
for k=3:(ds.nGammas + 2)
    indGammasCoarse(k) = indGammasCoarse(k-1) + ds.listCoarseGammas(k-2).nch;
end

% Second part: look for the close to touching parts using I and add those
% who belong to disc i. But just for the coarse discretization

disc_i = ds.I(1, 1);
disc_s = ds.I(1, 3);
k = 1;
while( (disc_i <= i || disc_s <= i) && k <= ds.nGammas )
    disc_i = ds.I(k, 1);
    disc_s = ds.I(k, 3);
    if(disc_i == i)
        % Meaning the first part of this gamma is on disc i
        endFlag = indGammasCoarse(k+1) + (ds.listCoarseGammas(k).nch)/2;
        flag( (indGammasCoarse(k+1)+1):endFlag, 1 ) = 1;
    end
    if( disc_s == i)
        % Meaning that the second part of this gamma is on disc i
        startFlag = indGammasCoarse(k+1) + (ds.listCoarseGammas(k).nch)/2 + 1;
        flag( startFlag:indGammasCoarse(k+2), 1 ) = 1;
    end
    % Update
    k = k+1;
end

end