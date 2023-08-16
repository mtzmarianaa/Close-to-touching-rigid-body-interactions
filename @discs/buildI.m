function [I, nGammas, pClose] = buildI(arg1, nDiscs)
% *buildI* builds the matrix storing information about the
% close-to-touching parts of a collection of non overlapping discs,
% $\Omega$. If there are $N_{\Omega}$ close-to-touching regions, then this
% matrix has $N_{\Omega}$ rows. At row $r$ this matrix has information
% regarding $\Gamma_k$ (close-to-touching part that involves disc i and
% disc s). Disc i's closest point to disc s is paramterized by
% $\theta_{ji}$. Disc s's closest point to disc i is parametrized by
% $\theta_{ks}$. Then row $r$ of I looks like:
%   $$ I_{r} = [i, \theta_{ji}, s, \theta_{ks}]. $$
%
% Syntax: [I, nGammas, pClose] = buildI(arg1, nDiscs)
%
% Input:
%   arg1 - either a discs object describing $\Omega$, collection of non overlapping discs or geom,
%             the description of the geometry of such discs (centers, radii)
%
% Optional input:
%   nDiscs - number of discs in the collection of non overlapping discs
%
% Output:
%   I - matrix describing the close-to-touching geometry
%   nGammas - number of close-to-touching regions
%   pClose - geometry of the close-to-touching regions (uniquely defined)
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

nGammas = 0;

% Determine if arg1 is a discs object or if its geom
if( isa(arg1, 'discs') && nargin < 2 )
    pClose = arg1.pClose;
    nDiscs = arg1.nDiscs;
else
    pClose = arg1;
end

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