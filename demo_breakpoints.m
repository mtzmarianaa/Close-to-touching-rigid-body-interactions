% Make a disc with given center and radius
center = [0 0];
radius = 1;
curve = @(t) disc(t, center=center, radius=radius);

% Manually specify panel breakpoints
breakpoints = linspace(0, 2*pi, 20);

% Make a chunker
chnkr = chunkerfuncbreakpoints(curve, breakpoints);

plot(chnkr, 'o-')
shg
