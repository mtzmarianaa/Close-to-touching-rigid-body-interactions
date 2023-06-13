% Verification of green's identity for discs (one disc first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


iseed = 166297;
rng(iseed,'twister');


fprintf("------n VERIFICATION OF GREEN'S IDENTITY ONE DISC-----\n\n\n ")

geom = [];
geom.ctrs = [0 ; 0];
geom.Rs = [1.0];
geom.nBreakPoints = 10;


ds = discs(geom);

figure(1)
plot(ds.chnkrs, '-o')
axis equal

figure(2)
plot(ds.chnkrs, '-o')
hold on
quiver(ds.chnkrs)
axis equal

%%%%% Place a point source f0 at x0 inside the disc, off-center but not too
% close to the boundary

x0 = 1/(2*sqrt(2))*[-1; -1];

%%%%% Evaluate u(x) = SL(x-x0) and DL(x - x0) at discretization points x

xOnSurface = reshape(ds.chnkrs.r, 2 , ds.chnkrs.k*ds.chnkrs.nch); % discretization points
nOnSurface = reshape(ds.chnkrs.n, 2 , ds.chnkrs.k*ds.chnkrs.nch);

% Generate off surface points randomly
tRand = 2*pi.*rand(1, 100);
normRand = 1.25 + (3)*rand(2, 100);
xOffSurface = normRand.*[sin(tRand); cos(tRand)];
xOnSurfaceTest = [sin(tRand); cos(tRand)];


% plot

figure(3)
plot(ds.chnkrs, 'b-o')
hold on
cq = [0 232/255 255/255];
quiver(ds.chnkrs, 'Color', cq)
cO = [147/255 155/255 255/255];
scatter(xOffSurface(1, :), xOffSurface(2, :), [], cO, "diamond")
cS = [155/255 0 255/255];
scatter(x0(1), x0(2), 75, cS, "o", 'filled')
axis equal

% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');


% Define u
u = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
gradSL = @(r) (-1)*r./(2*pi*vecnorm(r).*vecnorm(r));
pupn = @(n, x) dot(n, gradSL( bsxfun(@minus, x, x0) ));

% Plot u
figure(4)
scatter3(xOnSurface(1, :), xOnSurface(2, :), u(xOnSurface), 'Color', cq);
title('u on surface')

% Plot pupn
figure(5)
scatter3(xOnSurface(1, :), xOnSurface(2, :), pupn(nOnSurface, xOnSurface), 'Color', cO);
title('partial u partial n on surface')


% Use chunkerkerneval
DL_offSurface = chunkerkerneval(ds.chnkrs, DL_kern, u(xOnSurface), xOffSurface);
SL_offSurface = chunkerkerneval(ds.chnkrs, SL_kern, pupn(nOnSurface, xOnSurface), xOffSurface);
u_offSurface = u(xOffSurface);
DL_onSurface = chunkerkerneval(ds.chnkrs, DL_kern, u(xOnSurface), xOnSurfaceTest);
SL_onSurface = chunkerkerneval(ds.chnkrs, SL_kern, pupn(nOnSurface, xOnSurface), xOnSurfaceTest);
fprintf("%5.5e  : max u off surface: \n", max(abs(u_offSurface)))
u_onSurface = u(xOnSurface);
u_onSurfaceTest = u(xOnSurfaceTest);
fprintf("%5.5e  : max u on surface: \n", max(abs(u_onSurface)))
pupn_onSurface = pupn(nOnSurface, xOnSurface);
fprintf("%5.5e  : max pupn on surface: \n", max(abs(pupn_onSurface)))


% Plot DL, SL, u
figure(6)
plot(1:100, DL_offSurface, 'Color', cq, 'Marker', 'o')
hold on
plot(1:100, SL_offSurface, 'Color', cO, 'Marker', 'square')
plot(1:100, u_offSurface, 'Color', cS, 'Marker', '*')


err_OffSurface = DL_offSurface - SL_offSurface - u_offSurface';
err_OnSurface = (DL_onSurface + 0.5*u_onSurfaceTest') - SL_onSurface - u_onSurfaceTest';

figure(7)
plot(1:length(err_OffSurface), err_OffSurface, 'Color', cS);
title("Errors off surface")
fprintf("%5.5e  : max err off surface: \n", max(abs(err_OffSurface)))

figure(8)
plot(1:length(err_OnSurface), err_OnSurface, 'Color', cO);
title("Errors on surface")
fprintf("%5.5e  : max err on surface: \n", max(abs(err_OnSurface)))
