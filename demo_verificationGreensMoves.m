% Verification of green's identity for discs (one disc first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


iseed = 166297;
rng(iseed,'twister');
cq = [0 232/255 255/255];
cO = [147/255 155/255 255/255];
cS = [155/255 0 255/255];

fprintf("------n VERIFICATION OF GREEN'S IDENTITY ONE DISC-----\n\n\n ")

geom = [];
geom.ctrs = [0 ; 0];
geom.Rs = [1.0];
geom.nBreakPoints = 10;


ds = discs(geom);


%%%%% Get the discretization points of the boundary

xOnSurface = reshape(ds.chnkrs.r, 2 , ds.chnkrs.k*ds.chnkrs.nch);
nOnSurface = reshape(ds.chnkrs.n, 2 , ds.chnkrs.k*ds.chnkrs.nch);

% Generate off surface points randomly
tRand = 2*pi.*rand(1, 500);
normRand = 1.25 + (3)*rand(2, 500);
xOffSurface = normRand.*[sin(tRand); cos(tRand)];

% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');

%%%%% Place a point sources f0 at x0 inside the disc, off-center but not too
% close to the boundary

nx0 = linspace(2, 102, 101);

errs_OffSurface = zeros(size(nx0));
errs_OnSurface = zeros(size(nx0));

for i = 1:101
    % For each different norm of x0 verify Green's identity (we should see
    % a decrease in accurace as x0 moves closer to the boundary
    x0 = 1/(sqrt(2)*nx0(i))*[-1; -1];

    % Define u, gradSL, pupn
    u = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
    gradSL = @(r) (-1)*r./(2*pi*vecnorm(r).*vecnorm(r));
    pupn = @(n, x) dot(n, gradSL( bsxfun(@minus, x, x0) ));

    % Use chunkerkerneval
    DL_offSurface = chunkerkerneval(ds.chnkrs, DL_kern, u(xOnSurface), xOffSurface);
    SL_offSurface = chunkerkerneval(ds.chnkrs, SL_kern, pupn(nOnSurface, xOnSurface), xOffSurface);
    u_offSurface = u(xOffSurface);
    DLmat_onSurface = chunkermat(ds.chnkrs, DL_kern, []);
    SLmat_onSurface = chunkermat(ds.chnkrs, SL_kern, []);
    Imat = eye(size(DLmat_onSurface));
    u_onSurface = u(xOnSurface);
    pupn_onSurface = pupn(nOnSurface, xOnSurface);
    
    % Compute the l2 norm of the errors
    errs_OffSurface(i) = norm(DL_offSurface - SL_offSurface - u_offSurface');
    errs_OnSurface(i) = norm((0.5*Imat + DLmat_onSurface)*(u_onSurface') ...
        - SLmat_onSurface*(pupn_onSurface') - u_onSurface');

end


% Plot

figure(1)
plot(1./nx0, errs_OffSurface, 'Color', cS)
title("closeness to boundary vs l2 norm of errors off surface")


figure(2)
plot(1./nx0, errs_OnSurface, 'Color', cO)
title("closeness to boundary vs l2 norm of errors on surface")
