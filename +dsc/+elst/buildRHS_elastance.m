function [rhs_elastance, nu] = buildRHS_elastance(ds, listChunkers, nB, flagFunction, qk)
% *buildRHS_capacitance* builds the right hand side of the BIE for the
% capacitance problem. Depending on the input given it can build the RHS
% for the fine of the coarse discretizations.
%
% Syntax: rhs_capacitance = buildRHS_capacitance(ds, listChunkers, nB, flagFunction, uk)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   listChunkers - list of chunks of the geometry to use (either fine or
%                          coarse chunks on ds)
%   nB - stops for the chunks
%   flagFunction - function handle, which flag function to use - either flagnDisc or
%                          flagnDiscCoarse
%   qk - boundary data
%
% Output:
%   rhs_capacitance - rhs for the BIE (organized by blocks)
%   nu - $\nu$ in the elastance problem
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)


Ntot = nB(end);
nDiscs = ds.nDiscs;
chK = listChunkers(1).k;
nChunkers = length(listChunkers);

SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');

chnkrs = merge(listChunkers);


% Build nu (this is used for the RHS)
nuE = ones(nB(end), 1);
for i= 1:nDiscs
    % Use the flag
    flag = logical( flagFunction(i, ds) );
    % Convert that to points
    flag_points = repmat(flag, 1, chK);
    flag_points = flag_points';
    flag_points = logical( flag_points(:) );
    f = ones(nB(end), 1);
    f( ~flag_points) = 0; % Cancel out contribution from other discs
    % Get the perimeter of omegak
    per = chunkerintegral(chnkrs, f, []);
    nuE(flag_points) = qk(i)/per;
end

if(nargout>1)
    nu = nuE;
end

% Build the rhs matrix
rhs_chunkermat = zeros(Ntot);
opts = [];
opts2 = [];
opts2.adaptive_correction = true;

for k=1:nChunkers
    % Fill column by column
    chnkrk = listChunkers(k); % Current chunker (i.e. current gamma)
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nB(k) + 1;
    end_col = nB(k+1);
    for i=1:nChunkers
        start_row = nB(i) + 1;
        end_row = nB(i+1) ;
        chnkri = listChunkers(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i == k)
            % on boundary
            submat_rhs = chunkermat(chnkri, SL_kern, opts2);
        else
            % off boundary
            submat_rhs = chunkerkernevalmat(chnkri, SL_kern, targ, opts);
        end
        rhs_chunkermat(start_col:end_col, start_row:end_row) = submat_rhs;
    end
end

rhs_elastance = -rhs_chunkermat*nuE;

end
