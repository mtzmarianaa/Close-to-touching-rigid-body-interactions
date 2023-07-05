function [I, nGammas, pClose] = buildI(pClose, nDiscs)
% Builds the I matrix from pClose array of structs.
% At row r matrix I has information regarding gamma_r (close to touching
% part that involves disc i, disc s). Disc i's closest point to disc s is
% parametrized by theta_ji. Disc s's closest point to disc i is
% parametrized by theta_ks. Then row r of I looks like:
%   I(r) = [i, theta_ji, s, theta_ks]
%
% INPUT
%        pClose - information about the close to touching
%            regions (at least thetas, discClose, pRef given)
%                     pClose(i).data = (nClose x 3 matrix) at disc i the point
%                     closest to another disc in parameter
%                     space. The first column is the
%                     point in parameter (angle) space. The 
%                     second column is the index of
%                     the other disc which that point is
%                     closest to. The third column is the index
%                     of the point in the other disc that is
%                     closest to disc i.
%                     pClose(i).nClose = number of points at
%                     disc i considered to be in the close
%                     region
%                     pClose(i).thetasReg = angle of region
%                     considered to be close

nGammas = 0;

% For each disc we order the thetas and count the number of close to
% touching parts
for i=1:nDiscs
    [pClose(i).data, indSort] = sortrows( pClose(i).data );
    pClose(i).thetasReg = pClose(i).thetasReg(indSort);
    nClose = pClose(i).nClose;
    for j=1:nClose
        s = pClose(i).data(j, 2);
        if(s>i)
            % Meaning we haven't considered this gamma
            nGammas = nGammas + 1;
        end
    end
end

% Initialize I matrix
I = -1*ones(nGammas, 4);
k = 1;

% Now we build the matrix
for i=1:nDiscs
    nClose = pClose(i).nClose;
    for j=1:nClose
        s = pClose(i).data(j, 2);
        if(s>i)
            % Meaning we haven't added this gamma
            theta_ji = pClose(i).data(j, 1);
            ind_ks = pClose(i).data(j, 3);
            theta_ks = pClose(s).data( ind_ks, 1 );
            I(k, :) = [ i,  theta_ji, s, theta_ks];
            k = k+1;
        end
    end
end



end