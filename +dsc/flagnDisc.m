function [flag] = flagnDisc(i, ds)
% *flagnDisc* given a discs object outputs a sparse boolean vector such
% that it indicates which chunk on ds is on the i-th disc. Works only on
% the fine mesh. See flagnDiscCoarse if interested on the coarse mesh.
%
% Syntax: [flag] = flagnDisc(i, ds)
%
% Input:
%   i - index of discs we care about
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%
% Output:
%   flag - (ds.chnkrs.nch , 1) sparse boolean matrix, 0 if that chunk is
%            not on disc i, 1 if it is on disc i
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

assert( i <= ds.nDiscs, 'Not enough discs on ds for this request' )

flag = zeros(ds.chnkrs.nch, 1);

if( ~ds.infoClose )
    % Meaning that ds.listChnkrs is in disc order
    k = 1;
    nBk = zeros( i + 1, 1 );
    nBk(1) = 0;
    while(k < (i+1) )
        nBk(k+1) = nBk(k) + ds.listChnkrs(k).nch;
        k = k+1;
    end
    % Fill
    flag( nBk(i)+1:nBk(i+1), 1 ) = 1;
else
    % Meaning that ds.listChnkrs has information from close to touching
    % parts and it's not in disc order
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
    % Second part: look for the close to touching parts using I and add
    % those who belong to disc i
    disc_i = ds.I(1, 1);
    disc_s = ds.I(1, 3);
    k = 1;
    while( (disc_i <= i || disc_s <= i) && k <= ds.nGammas )
        disc_i = ds.I(k, 1);
        disc_s = ds.I(k, 3);
        if(disc_i == i)
            % Meaning the first part of this gamma is on disc i
            endFlag = ds.indGammas(k+1) + (ds.listGammas(k).nch)/2;
            flag( (ds.indGammas(k+1)+1):endFlag, 1 ) = 1;
        end
        if( disc_s == i)
            % Meaning that the second part of this gamma is on disc i
            startFlag = ds.indGammas(k+1) + (ds.listGammas(k).nch)/2 + 1;
            flag( startFlag:ds.indGammas(k+2), 1 ) = 1;
        end
        % Update
        k = k+1;
    end
end

end