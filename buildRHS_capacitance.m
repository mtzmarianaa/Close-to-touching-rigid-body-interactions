function rhs_capacitance = buildRHS_capacitance(ds, listChunkers, nB, flagFunction, uk)
% Build the rhs to solve the capacitance problem. Depending on the arguments
% given it can build the rhs for the fine or the coarse discretizations

Ntot = nB(end);
nDiscs = ds.nDiscs;
chK = listChunkers(1).k;

chnkrs = merge(listChunkers);

% Fill the RHS
rhs_capacitance = zeros(Ntot, 1);
for i = 1:nDiscs 
    % Use the flag
    flag = logical( flagFunction(i, ds) );
    % Convert that to points
    flag_points = repmat(flag, 1, chK);
    flag_points = flag_points';
    flag_points = logical( flag_points(:) );
    % Get the points
    xOnSurface = reshape(chnkrs.r(:, :, flag ), 2, []);
    u_toUse = uk{i}(xOnSurface);
    rhs_capacitance(flag_points) = u_toUse';
end

end