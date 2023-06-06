% Test the discs class

clear all
close all

geom = [];
geom.ctrs = [-1 0 3 7; -6 0 2 4];
geom.Rs = [1, 1.3, 0.4, 1.2];

ds = discs(geom);

figure(1)
plot(ds.chnkrs, '-o')
axis equal

figure(2)
plot(ds.chnkrs, '-o')
hold on
quiver(ds.chnkrs)
axis equal

