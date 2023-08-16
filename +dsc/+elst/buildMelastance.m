function M = buildMelastance(ds, listChunkers, nB, flagFunction)
% *buildMelastance* builds the M matrix necessary for the elastanec
% problem. Depending on the input given it can build either the M for the
% fine discretization or the M for the coarse mesh. This is a block
% diagonal matrix.
%
% Syntax: M = buildMelastance(ds, listChunkers, nB, flagFunction)
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
%
% Output:
%   M - (ds.chnkrs.npt, ds.chnkrs.npt) M matrix
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

Ntot = nB(end);
nDiscs = ds.nDiscs;
chK = listChunkers(1).k;
nChunkers = length(listChunkers);

M = zeros(Ntot);
all_weights = zeros(1, Ntot);

% Fill in the all weights vector
for i=1:nChunkers
    w = weights(listChunkers(i));
    all_weights( (nB(i)+1):nB(i+1) ) = reshape(w, listChunkers(i).npt, [])';
end

% Fill in the matrix
for i=1:nChunkers
    start_index = nB(i) + 1;
    end_index = nB(i+1);
    % Now filter according to discs
    for k=1:nDiscs
        flag = logical( flagFunction(k, ds) );
        % Convert that to points
        flag_points = repmat(flag, 1, chK);
        flag_points = flag_points';
        flag_points = logical( flag_points(:) );
        flag_points_here = false( 1, Ntot);
        flag_points_here(start_index:end_index) = flag_points(start_index:end_index);
        % Fill the matrix
        weights_here = zeros( 1, Ntot );
        weights_here(flag_points') = all_weights(flag_points');
        weights_here = repmat(weights_here, sum(flag_points(start_index:end_index)), 1);
        M( flag_points_here, :) = weights_here;
    end
end

end