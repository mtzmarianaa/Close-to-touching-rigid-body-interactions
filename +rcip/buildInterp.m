function matInterpolant = buildInterp(listPrecomputedR, uDist)
% Given a list of precomputed R matrices build a nxnx16 matrix with
% coefficients for interpolation

kDist = size(listPrecomputedR{1}, 3); % Get the degree of the Legendre polynomial wanted
nCh = size(listPrecomputedR, 2); % Number of chunks (distance chunks)
nR = size(listPrecomputedR{1}, 1); 

if(nargin < 2)
    [~ , ~ , uDist , ~] = lege.exps(kDist);
end

matInterpolant = cell(1, nCh); % Initialize cell array
cellItem = zeros(size(listPrecomputedR{1})); % Initialize item on cell array

for ch=1:nCh
    % For each chunk select a group of matrices
    Rch = listPrecomputedR{ch};
    for row=1:nR
        for col=1:nR
            q = reshape(Rch(row, col, :), kDist, 1);
            % Compute matrix multiplication
            coefs = uDist*q;
            % Save
            cellItem(row, col, :) = coefs;
        end
    end
    matInterpolant{ch} = cellItem;
end

end