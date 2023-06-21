% Verification of green's identity for discs (one disc first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all



fprintf("------n VERIFICATION OF GREEN'S IDENTITY ONE DISC-----\n\n\n ")

geom = [];
geom.ctrs = [0  ; 0 ];
geom.Rs = [2.0 ];
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

%x0 = 1/(2*sqrt(2))*[1; 1];

x0 = [0;0];

%%%%% Evaluate u(x) = SL(x-x0) and DL(x - x0) at discretization points x

xOnSurface = reshape(ds.chnkrs.r, 2 , ds.chnkrs.k*ds.chnkrs.nch); % discretization points
nOnSurface = reshape(ds.chnkrs.n, 2 , ds.chnkrs.k*ds.chnkrs.nch);

% Generate off surface points randomly
rng(166297);
tRand = 2*pi.*rand(1, 100);
normRand = 1.25 + (3)*rand(1, 100);
normRand0 = 1.1 + 0.5*rand(1, 100);
rng(257);
tRand0 = 2*pi.*rand(1, 100);
xOffSurface = [normRand0.*sin(tRand) normRand0.*sin(tRand0) + 3; normRand.*cos(tRand) normRand.*cos(tRand0)];
nOffSurface = size(xOffSurface, 2);

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


% Define u, gradSL, pupn
u = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
gradSL = @(r) (-1)*r./(2*pi*vecnorm(r).*vecnorm(r));
pupn = @(n, x) dot(n, gradSL( bsxfun(@minus, x, x0) ));

% Plot u
figure(4)
scatter3(xOnSurface(1, :), xOnSurface(2, :), u(xOnSurface), 'Color', cq);
title('u on surface')
xlim([-2, 5])
ylim([-2, 5])

% Plot pupn
figure(5)
scatter3(xOnSurface(1, :), xOnSurface(2, :), pupn(nOnSurface, xOnSurface), 'Color', cO);
title('partial u partial n on surface')
xlim([-2, 5])
ylim([-2, 5])

% Use chunkerkerneval
DL_offSurface = chunkerkerneval(ds.chnkrs, DL_kern, u(xOnSurface), xOffSurface);
SL_offSurface = chunkerkerneval(ds.chnkrs, SL_kern, pupn(nOnSurface, xOnSurface), xOffSurface);
u_offSurface = u(xOffSurface);
DLmat_onSurface = chunkermat(ds.chnkrs, DL_kern, []);
SLmat_onSurface = chunkermat(ds.chnkrs, SL_kern, []);
Imat = eye(size(DLmat_onSurface));
fprintf("%5.5e  : max u off surface: \n", max(abs(u_offSurface)))
u_onSurface = u(xOnSurface);
fprintf("%5.5e  : max u on surface: \n", max(abs(u_onSurface)))
pupn_onSurface = pupn(nOnSurface, xOnSurface);
fprintf("%5.5e  : max pupn on surface: \n", max(abs(pupn_onSurface)))


% Plot DL, SL, u
figure(6)
plot(1:nOffSurface, DL_offSurface, 'Color', cq, 'Marker', 'o')
hold on
plot(1:nOffSurface, SL_offSurface, 'Color', cO, 'Marker', 'square')
plot(1:nOffSurface, u_offSurface, 'Color', cS, 'Marker', '*')


err_OffSurface = DL_offSurface - SL_offSurface - u_offSurface';
err_OnSurface = (0.5*Imat + DLmat_onSurface)*(u_onSurface') - SLmat_onSurface*(pupn_onSurface') - u_onSurface';
fprintf("%5.5e : max first term\n", max(abs((0.5*Imat + DLmat_onSurface)*(u_onSurface') )) )
fprintf("%5.5e : max second term\n", max(abs( SLmat_onSurface*(pupn_onSurface') )) )
fprintf("%5.5e : max third term\n", max(abs( u_onSurface' )) )

figure(7)
plot(1:nOffSurface, err_OffSurface, 'Color', cS);
title("Errors off surface")
fprintf("%5.5e  : max err off surface: \n", max(abs(err_OffSurface)))
fprintf("%5.5e  : norm err off surface: \n", norm(err_OffSurface) )

figure(8)
plot(1:length(err_OnSurface), err_OnSurface, 'Color', cO);
title("Errors on surface")
fprintf("%5.5e  : max err on surface: \n", max(abs(err_OnSurface)))
fprintf("%5.5e  : norm err on surface: \n", norm(err_OnSurface) )










