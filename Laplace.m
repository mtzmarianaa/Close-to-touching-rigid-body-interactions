% Laplace problem (1 disc first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

cq = [0 232/255 255/255];
cO = [147/255 155/255 255/255];
cS = [155/255 0 255/255];
%%%%%%%%%

% For a point charge
x0 = 0.25/(2*sqrt(2))*[1; 1];

u1 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
u2 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);

uk = {u1};
ctrs = [0 ;0 ];
Rs = [1.0];
n = length(uk);
nBreakPoints = [10];


% Define points on surface
geom = [];
geom.ctrs = ctrs;
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
ds = discs(geom);
nk = (ds.nBreakPoints - 1).*16; % Number of discretization points on each disk
nk = [0; nk]; % Useful for filling in the matrix K
Ntot = sum(nk);

% Plot the circles
figure(1)
plot(ds.chnkrs, '-o', 'Color', cq)
hold on
quiver(ds.chnkrs, 'Color', cO)
scatter(x0(1), x0(2), 'Color', cS)
axis equal


% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');
DLplusSL = @(s,t) DL_kern(s,t) + SL_kern(s,t);

% Initialize the matrix
K_lap = zeros(Ntot);
opts = [];

% Proceed to fill in the matrix 
for k=1:n
    % fill column by column
    nCol = nk(k + 1); % Number of columns for the submatrices
    chnkrk = ds.listChnkrs(k);
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nk(k) + 1;
    end_col = nk(k) + nCol;
    for i=1:n
        nRow = nk(i + 1); % Number of rows for the submatrix
        start_row = nk(i) + 1;
        end_row = nk(i) + nRow;
        chnkri = ds.listChnkrs(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i == k)
            % on boundary
            submat = chunkermat(chnkri, DLplusSL) + 0.5*eye(nRow);
        else
            % off boundary
            submat = chunkerkernevalmat(chnkri, DLplusSL, targ, opts);
        end
        % Add this submatrix to K
        K_lap(start_row:end_row, start_col:end_col) = submat;
    end
end

% Fill the RHS

rhs_lap = zeros(Ntot, 1);

for i = 1:n
    start_row = nk(i) + 1;
    end_row = nk(i) + nk(i + 1);
    chnkri = ds.listChnkrs(i); 
    xOnSurface = reshape(chnkri.r, 2, chnkri.k*chnkri.nch);
    u_toUse = uk{i}(xOnSurface);
    rhs_lap(start_row:end_row) = u_toUse;
end


% Solve for sigma, the unknown density

s = tic();
sigma_lap = K_lap\rhs_lap;
t1 = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with Matlab's backslash (LAPLACE)\n", t1);


%%%%% Plot u off surface
% Find points off surface
rmin = min(ds.chnkrs ) - 5*[1;1]; 
rmax = max(ds.chnkrs ) + 5*[1;1];
nplot = 250;
hx = (rmax(1)-rmin(1))/nplot;
hy = (rmax(2)-rmin(2))/nplot;
xtarg = linspace(rmin(1)+hx/2,rmax(1)-hx/2,nplot); 
ytarg = linspace(rmin(2)+hy/2,rmax(2)-hy/2,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

s = tic; 
in = chunkerinterior(ds.chnkrs , targets ); 
t3 = toc(s);
fprintf('%5.2e s : time to find points off surface\n',t3)

% Evaluate points off surface
s = tic;
Dsol = chunkerkerneval(ds.chnkrs , DLplusSL,  sigma_lap , targets(:,~in)); 
t4 = toc(s);
fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t4);

% Plot off surface
figure(7)
zztarg = nan(size(xxtarg));
zztarg(~in) = Dsol;
h=surf(xxtarg,yytarg,zztarg);
set(h,'EdgeColor','none')
title("Solution for exterior Laplace")
colorbar