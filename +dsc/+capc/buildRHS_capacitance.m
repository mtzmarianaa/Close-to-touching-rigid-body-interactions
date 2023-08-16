function rhs_capacitance = buildRHS_capacitance(ds, listChunkers, nB, flagFunction, uk)
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
%   uk - boundary data
%
% Output:
%   rhs_capacitance - rhs for the BIE (organized by blocks)
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

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