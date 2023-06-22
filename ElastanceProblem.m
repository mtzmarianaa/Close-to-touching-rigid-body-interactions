% Elastance problem (2 discs first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%close all

% Define the charges on each disc
q1 = -4.162781873389288733733337721787393093109130859375;
q2 = 1.2219951298941718054180682884179987013339996337890625;

qs = [q1; q2];


% Define the geometry of the circles

ctrs = [0 5;0 0]; % Centers of the circles
Rs = [1.5; 1.25]; % Radi of the circles
n = length(qs);
nBreakPoints = [10; 10];


% Define points on surface
geom = [];
geom.ctrs = ctrs;
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
ds = discs(geom);
Ntot = ds.chnkrs.npt;
nB = (ds.nBreakPoints - 1).*16; % Number of discretization points on each disk
nB = [0 nB]; % Useful for filling in the matrix K
nk = zeros(size(nB));
for i=2:(n+1)
    nk(i) = nk(i-1) + nB(i);
end

% Plot the circles
figure()
cq = [0 232/255 255/255];
cO = [147/255 155/255 255/255];
cS = [155/255 0 255/255];
plot(ds.chnkrs, '-o', 'Color', cq)
hold on
quiver(ds.chnkrs, 'Color', cO)
axis equal
title("Elastance Problem - Circles and Normals")

% Plot the charges qk
figure()
for i=1:n
    [X,Y,Z] = ellipsoid( ctrs(1, i) , ctrs(2, i) , qs(i) , Rs(i) , Rs(i) , 0 );
    s = surf(X,Y,Z,'FaceAlpha',0.5);
    s.EdgeColor = 'none';
    hold on
end
colorbar
title("Elastance Problem - Total charge on each disc, given")
hold off


% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');
DLplusSL = @(s,t) DL_kern(s,t) + SL_kern(s,t);
ID_kern = @(s,t)  speye(size(t.r, 2)*size(t.r, 3), size(s.r, 2) *size(s.r, 3));

% Define the points on surface
xOnSurface = reshape(ds.chnkrs.r, 2, ds.chnkrs.k*ds.chnkrs.nch);

% Build nu (this is used for the RHS)
nu = ones(ds.chnkrs.npt, 1);
for k= 1:n
    chnkrk = ds.listChnkrs(k);
    nPointsk = chnkrk.npt;
    start_row = nk(k) + 1;
    end_row = nk(k + 1);
    f = ones(nPointsk, 1);
    % Get the perimeter of omegak
    per = chunkerintegral(chnkrk, f, []);
    nu(start_row:end_row) = qs(k)/per;
end

% Build the RHS matrix, the I/2 + DL matrix and the M matrix
IDplusDLMat = zeros(Ntot);
opts = [];

for k=1:n
    % Fill column by column
    nCol = nk(k+1) - nk(k); % Number of columns for the submatrices
    chnkrk = ds.listChnkrs(k); % Current circle
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nk(k) + 1;
    end_col = nk(k+1);
    for i=1:n
        nRow = nk(i + 1) - nk(i); % Number of rows for the submatrix
        start_row = nk(i) + 1;
        end_row = nk(i+1) ;
        chnkri = ds.listChnkrs(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i == k)
            % on boundary
            submat = chunkermat(chnkri, DL_kern) + 0.5*eye(nRow);
        else
            % off boundary
            submat = chunkerkernevalmat(chnkri, DL_kern, targ, opts);
        end
        % Add this submatrix to K
        IDplusDLMat(start_col:end_col, start_row:end_row) = submat;
    end
end

% Build the block diagonal matrix M
M = zeros(Ntot);
for i=1:n
    start_index = nk(i) + 1;
    end_index = nk(i+1);
    chnkri = ds.listChnkrs(i); 
    w = weights(chnkri);
    w = reshape(w, chnkri.npt, [])';
    submat = repmat(w, chnkri.npt, 1); % Block matrix with all its rows being the weights
    M(start_index:end_index, start_index:end_index) = submat;
end

% Complete building the rhs 
rhs_chunkermat = chunkermat(ds.chnkrs, SL_kern,  []);
rhs = -rhs_chunkermat*nu;


% Build the system
%IDplusDLMat = chunkermat(ds.chnkrs, DL_kern, opts);
%IDplusDLMat = IDplusDLMat + 0.5*speye(size(IDplusDLMat));
K = IDplusDLMat + M;

% Solve for sigma, unknown density
s = tic();
sigma = K\rhs;
t = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with Matlab's backslash\n", t);


% Compute u for on surface points

u_onSurfaceElastance = -M*sigma;

% Plot u on surface
figure()
scatter3(xOnSurface(1, :), xOnSurface(2, :), u_onSurfaceElastance, 'Color', cq);
title('Elastance Problem - u on surface')
xlim([-2, 8])
ylim([-2, 8])


% Plot u off surface


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
Dsol_elastance = chunkerkerneval(ds.chnkrs , DL_kern,  sigma , targets(:,~in));
Dsol_elastance = Dsol_elastance + chunkerkerneval(ds.chnkrs, SL_kern, nu, targets(:, ~in));
t4 = toc(s);
fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t4);

% Plot off surface
figure()
zztarg = nan(size(xxtarg));
zztarg(~in) = Dsol_elastance;
h=surf(xxtarg,yytarg,zztarg);
set(h,'EdgeColor','none')
title("Elastance Problem - Potential off surface")
colorbar



% If capacitance problem was run before
difSolOff = Dsol_elastance - Dsol_capacitance;
figure()
zztarg = nan(size(xxtarg));
zztarg(~in) = difSolOff;
h=surf(xxtarg,yytarg,zztarg);
set(h,'EdgeColor','none')
title("Diference in potentials Capacitance and Elastance")
colorbar
fprintf('%5.2e s : Max error between computed potentials\n',max(abs(difSolOff)));
fprintf('%5.2e s : l2 error between computed potentials\n', norm(difSolOff));

