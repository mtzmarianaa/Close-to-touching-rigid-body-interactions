function ds = splitCloseRegions(ds, opts)
% Given a discs object, splits those chunks that are considered to be in
% the close to touching region

nDiscs = ds.nDiscs;

if(nargin<2)
    opts = [];
end

[x, w, u] = lege.exps( ds.chnkrs.k ); % precompute

% For each disc check if there are close to touching regions

for i=1:nDiscs
    chClosei = ds.indCloseChunk(i).ind; % We know that they are in order
    if( ~isempty(chClosei) )
        % Meaning we have close to touching regions
        nClosei = length(chClosei);
        chnkri = ds.listChnkrs(i); % Current chunker object we want to split
        offset = 0:(nClosei - 1);  % Offset of the new chunks vs the original ones
        newClose = [chClosei + offset'; (chClosei + offset' + 1)];
        for k=1:nClosei
            % Split each close chunk into two close chunks
            ich = chClosei(k) + offset(k);
            chnkri = split(chnkri, ich, opts, x, w, u); % Split current chunker
        end
        % Update
        ds.indCloseChunk(i).ind = sort(newClose); % Assume everything contained in the original close region is close
        ds.listChnkrs(i) = chnkri;
        ds.nBreakPoints = ds.nBreakPoints + nClosei;
    end
end

% Update the big chunker object
ds.chnkrs = merge(ds.listChnkrs);

end