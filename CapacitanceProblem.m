% Capacitance problem (2 discs first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% For a point charge
x0 = 1/(2*sqrt(2))*[1; 1];

u1 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
u2 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);


uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 3;0 0]; % Centers of the circles
Rs = [1.0; 1.0]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];


% Define points on surface
geom = [];
geom.ctrs = ctrs;
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
ds = discs(geom);
nk = (nBreakPoints - 1).*16; % Number of discretization points on each disk
nk = [0; nk]; % Useful for filling in the matrix K
Ntot = sum(nk);

% Plot the circles
figure(1)
cq = [0 232/255 255/255];
cO = [147/255 155/255 255/255];
cS = [155/255 0 255/255];
plot(ds.chnkrs, '-o', 'Color', cq)
hold on
quiver(ds.chnkrs, 'Color', cO)
axis equal

% Plot the uks on surface
figure(2)
for i =1:n
    chnkri = ds.listChnkrs(i); 
    xOnSurface = reshape(chnkri.r, 2, chnkri.k*chnkri.nch);
    u_toUse = uk{i}(xOnSurface);
    scatter3(xOnSurface(1, :), xOnSurface(2, :), u_toUse, [], cO)
    hold on
end
title("uk on surface")
hold off


% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');
DLplusSL = @(s,t) DL_kern(s,t) + SL_kern(s,t);

% Initialize the matrix
K = zeros(Ntot);
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
        K(start_row:end_row, start_col:end_col) = submat;
    end
end

% Fill the RHS

rhs = zeros(Ntot, 1);

for i = 1:n
    start_row = nk(i) + 1;
    end_row = nk(i) + nk(i + 1);
    chnkri = ds.listChnkrs(i); 
    xOnSurface = reshape(chnkri.r, 2, chnkri.k*chnkri.nch);
    u_toUse = uk{i}(xOnSurface);
    rhs(start_row:end_row) = u_toUse;
end


% Solve for sigma, the unknown density

s = tic();
sigma_B = K\rhs;
t1 = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with Matlab's backslash\n", t1);

s = tic();
sigma_G = gmres(K, rhs, [], 1e-14, 100);
t2 = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with GMRES\n", t2);


% Plot the sigma_ks on surface
figure(3)
allXonSurface = reshape(ds.chnkrs.r, 2, ds.chnkrs.k*ds.chnkrs.nch);
scatter3(allXonSurface(1, :), allXonSurface(2, :), sigma_B, [], cq)
title("sigma on surface, Matlab's backslash")


% Plot the sigma_ks on surface
figure(4)
scatter3(allXonSurface(1, :), allXonSurface(2, :), sigma_G, [], cq)
title("sigma on surface, GMRES")



% Compute the boundary integral of sigma over the discs to find q
q_B = zeros(n, 1);
q_G = zeros(n, 1);
opts.usesmooth = false;

for i=1:n
    chnkri = ds.listChnkrs(i);
    start_row = nk(i) + 1;
    end_row = nk(i) + nk(i + 1);
    sigma_Bi = sigma_B(start_row:end_row);
    q_B(i) = chunkerintegral(chnkri, sigma_Bi, opts);
    fprintf("%5.5e : charge at disk %1.0e with backslash\n", q_B(i), i);
    sigma_Gi = sigma_G(start_row:end_row);
    q_G(i) = chunkerintegral(chnkri, sigma_Gi, opts);
    fprintf("%5.5e : charge at disk %1.0e with GMRES\n\n", q_G(i), i);
end



% Plot the charges qk
figure(5)
for i=1:n
    [X,Y,Z] = ellipsoid( ctrs(1, i) , ctrs(2, i) , q_B(i) , Rs(i) , Rs(i) , 0 );
    s = surf(X,Y,Z,'FaceAlpha',0.5);
    s.EdgeColor = 'none';
    hold on
end
colorbar
title("Total charge on each disc, Matlab's backslash")
hold off



% Plot the charges qk
figure(6)
for i=1:n
    [X,Y,Z] = ellipsoid( ctrs(1, i) , ctrs(2, i) , q_G(i) , Rs(i) , Rs(i) , 0 );
    s = surf(X,Y,Z,'FaceAlpha',0.5);
    s.EdgeColor = 'none';
    hold on
end
colorbar
title("Total charge on each disc, GMRES")
hold off

