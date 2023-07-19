function M = buildMelastance(ds, listChunkers, nB, flagFunction)
% Builds the M function necessary for the elastance problem. Depending on
% the arguments given it can build either the M for the fine discretization
% or the M for the coarse discreatization.
% Build the block diagonal matrix M (note that here we don't care about the
% order of the discs, just the order of the chunker objects - order of
% either discs (if no close information) or the order of gammas (if close
% information)

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