breakpoints1 = linspace(0, 2*pi, 10);
breakpoints2 = linspace(0, 2*pi, 20);

chnkr1 = chunkerfuncbreakpoints(@(t) disc(t, center=[0 0], radius=1),   breakpoints1);
chnkr2 = chunkerfuncbreakpoints(@(t) disc(t, center=[3 2], radius=0.5), breakpoints2);
chnkrs = merge([chnkr1 chnkr2]);

plot(chnkrs, 'o-')
shg
