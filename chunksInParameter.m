function chs = chunksInParameter(chnkr, thetas, x)
% Given a chunker object get all the chunks in parameter space
%
% Input:
%        chnkr - the chunker object
%        thetas - breakpoints in the parameter space

k = chnkr.k;
nch = chnkr.nch; % Number of chunks
chs = zeros(k, nch);


if( nargin < 3 )
    x = lege.exps(chnkr.k);
end


if( nargin < 2)
    thetas = linspace(0, 2*pi, nch);
end

ths = [thetas thetas(1)]; % For the loop

% Do this for all the chunks
for i = 1:nch
    tha = ths(i);
    thb = ths(i+1);
    chs(:, i) = chunkInParameter(chnkr, tha, thb, x);
end

end