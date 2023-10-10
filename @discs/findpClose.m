function pClose = findpClose(ds, maxDist, thetaReg)
% *findpClose*  given a discs object compute the closest points and store
% it in a struct called pClose.
%
% Syntax: pClose = findpClose(ds)
%              pClose = findpClose(ds, maxDist)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs
%
% Optional input:
%   maxDist - max distance between discs to be considered as "close" (0.1)
%
% Output:
%   pClose - list of structs containing the close-to-touching parts of the discs
%                 pClose(i).data = [thetaMin, neighbor index, number of close to touching]
%                 pClose(i).nClose = number of discs close to disc i
%                 pClose(i).thetasReg = +- theta for the close-to-region
%                                                    part
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

if(nargin<2)
    maxDist = 0.1;
end

if(nargin<3)
    thetaReg = pi/6;
end

pClose = [];
nDiscs = ds.nDiscs; % Get number of discs in this collection of discs
currentN = ones(nDiscs, 1); % Current index of the neighbors on each disc
fullInfo = zeros(6, 3, nDiscs);

% Iterate through the discs
for i=1:nDiscs
    k = 1;
    neis = zeros(6, 1); % Maximum number of neighbors possible for a disc
    % First get the index of the neighbors
    for j=1:nDiscs
        if(j~=i)
            if( norm(ds.ctrs(:, i)- ds.ctrs(:, j)) - ds.Rs(i) - ds.Rs(j) <= maxDist )
                neis(k) = j;
                if(k == 6)
                    break;
                end
                k = k+1;
            end
        end
    end
    % Remove those rows with zeros in the neis vector
    neis = neis(any(neis,2), :);
    pClose(i).neis = neis;
    pClose(i).nClose = length(neis);
    pClose(i).data = zeros(pClose(i).nClose, 3);
    pClose(i).thetasReg = thetaReg;

end


for i =1:nDiscs
    % Then get the data for the close to touching interactions,
    % pClose(i).data but we just iterate through the neighbors
    for j=1:pClose(i).nClose
        s = pClose(i).neis(j); % Current neighbor we are considering
        nei_i = currentN(i);
        nei_s = currentN(s);
        if( i<s )
            % Otherwise we've already considered this close-to-touching
            % part
            p_i = ds.Rs(i)*(ds.ctrs(:, s) - ds.ctrs(:, i))/norm(ds.ctrs(:, s) - ds.ctrs(:, i)) + ds.ctrs(:, i);
            p_s = ds.Rs(s)*(ds.ctrs(:, i) - ds.ctrs(:, s))/norm(ds.ctrs(:, i) - ds.ctrs(:, s)) + ds.ctrs(:, s);
            % We need the close-to-touching point in parameter space
            % For disc i
            if(norm(ds.ctrs(:, i)) > 1e-12)
                % The center is not zero
                th_i = (ds.ctrs(:, s) - ds.ctrs(:, i))'*ds.ctrs(:, i);
                th_i = th_i/norm(th_i);
                th_i = acos( th_i );
            else
                % The center is zero
                th_i = acos( p_i(1)/norm(p_i) );
            end
            % For disc s
            if(norm(ds.ctrs(:, s)) > 1e-12)
                % The center is not zero
                th_s = (ds.ctrs(:, i) - ds.ctrs(:, s))'*ds.ctrs(:, s);
                th_s = th_s/norm(th_s);
                th_s = acos( th_s );
            else
                % The center is zero
                th_s = acos( p_s(1)/norm(p_s) );
            end
            % Save
            pClose(i).data(nei_i, :) = [th_i s nei_s];
            pClose(s).data(nei_s, :) = [th_s i nei_i];
            currentN(i) = currentN(i) + 1;
            currentN(s) = currentN(s) + 1;
        end
    end
end

% Remove the unnecessary neis attribute from pClose
pClose = rmfield(pClose,'neis');

end