% Test the discs class

clear all
close all



fprintf("------n EVALUATING INTEGRAL OF SMOOTH INTEGRANDS-----\n\n\n ")

%%%%%%%%% HUGE, NON ADAPTIVE

geom = [];
geom.ctrs = 25*[-1 0 3 7; -6 0 2 4];
geom.Rs = 25*[1, 1.3, 0.4, 1.2];
geom.nBreakPoints = 3;

ds = discs(geom);

figure(1)
plot(ds.chnkrs, '-o')
axis equal

figure(2)
plot(ds.chnkrs, '-o')
hold on
quiver(ds.chnkrs)
axis equal


% We try to use chunkie to evaluate smooth integrands from values given at
% the surface discretization nodes
points = reshape(ds.chnkrs.r, 2 , ds.chnkrs.k*ds.chnkrs.nch);

f2 = @(x) f1_smooth(x);
f2vals = f2(points);


% Just for sanity check we plot these
figure(3)
scatter3(points(1, :), points(2, :), f2vals)

opts = [];
opts.usesmooth = true;

int1 = chunkerintegral(ds.chnkrs, f2vals, opts);
true_int1 = 2*pi*(sum(geom.Rs));
fprintf("-----HUGE CIRCLES, NON ADAPTIVE 3 BREAKPOINTS, AUTOMATIC PREF ---\n\n\n")
fprintf("\n\n\n%5.15e  : computed integral.\n\n", int1)
fprintf("%5.5e  : absolute error.\n\n", abs(int1 - true_int1))
fprintf("%5.5e  : relative error.\n\n", abs(int1 - true_int1)/int1 )




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% HUGE, ADAPTIVE

geom = [];
geom.ctrs = 25*[-1 0 3 7; -6 0 2 4];
geom.Rs = 25*[1, 1.3, 0.4, 1.2];
geom.nBreakPoints = 10;

ds2 = discs(geom);

figure(4)
plot(ds2.chnkrs, '-o')
axis equal

figure(5)
plot(ds2.chnkrs, '-o')
hold on
quiver(ds2.chnkrs)
axis equal


% We try to use chunkie to evaluate smooth integrands from values given at
% the surface discretization nodes
points = reshape(ds2.chnkrs.r, 2 , ds2.chnkrs.k*ds2.chnkrs.nch);

f2 = @(x) f1_smooth(x);
f2vals = f2(points);


% Just for sanity check we plot these
figure(6)
scatter3(points(1, :), points(2, :), f2vals)

opts = [];
opts.usesmooth = true;

int1 = chunkerintegral(ds2.chnkrs, f2vals, opts);
true_int1 = 2*pi*(sum(geom.Rs));
fprintf("-----HUGE CIRCLES, NON ADAPTIVE 10 BREAKPOINTS, AUTOMATIC PREF ---\n\n\n")
fprintf("\n\n\n%5.15e  : computed integral.\n\n", int1)
fprintf("%5.5e  : absolute error.\n\n", abs(int1 - true_int1))
fprintf("%5.5e  : relative error.\n\n", abs(int1 - true_int1)/int1 )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% HUGE, ADAPTIVE, ANOTHER FUNCTION  

geom = [];
geom.ctrs = [-1 0 3 7; -6 0 -2 4];
geom.Rs = [1, 1.3, 0.8, 1.2];
geom.nBreakPoints = 20;


ds3 = discs(geom);

figure(7)
plot(ds3.chnkrs, '-o')
axis equal

figure(8)
plot(ds3.chnkrs, '-o')
hold on
quiver(ds3.chnkrs)
axis equal


% We try to use chunkie to evaluate smooth integrands from values given at
% the surface discretization nodes
points = reshape(ds3.chnkrs.r, 2 , ds3.chnkrs.k*ds3.chnkrs.nch);

f2 = @(x) f2_smooth(x);
f2vals = f2(points);


% Just for sanity check we plot these
figure(9)
scatter3(points(1, :), points(2, :), f2vals)

opts = [];
opts.usesmooth = true;

int1 = chunkerintegral(ds3.chnkrs, f2vals, opts);
fprintf("-----HUGE CIRCLES, NON ADAPTIVE 10 BREAKPOINTS\n\n\n")
fprintf("\n\n\n%5.15e  : computed integral.\n\n", int1)


