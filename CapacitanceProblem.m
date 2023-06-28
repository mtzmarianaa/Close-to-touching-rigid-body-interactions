function [q, sigma, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, plt, outopt)
% Given the description of the geometry of nCircles solve the capacitance
% problem. 
% IN: ds : a discs object with the given geometry of the discs
%        uk:  u defined at the boundary
%        plot: boolean, wheather to plot or not
% OUT:  q: charges on the n given discs
%           sigma: on bundary density

if( nargin < 3 )
    plt = false;
end
if( nargin < 4 )
    outopt = [];
    outopt.fact = 3;
    outopt.nplot = 250;
end

if ~isfield(outopt, 'fact')
    outopt.fact = 3;
end


if ~isfield(outopt, 'nplot')
    outopt.nplot = 250;
end

n = length(uk);
ctrs = ds.ctrs;
Rs = ds.Rs;
nBreakPoints = ds.nBreakPoints;


Ntot = ds.chnkrs.npt;
nB = (ds.nBreakPoints - 1).*16; % Number of discretization points on each disk
nB = [0 nB]; % Useful for filling in the matrix K
nk = zeros(size(nB));
for i=2:(n+1)
    nk(i) = nk(i-1) + nB(i);
end



if( plt )
    % Plot the circles
    figure()
    cq = [0 232/255 255/255];
    cO = [147/255 155/255 255/255];
    cS = [155/255 0 255/255];
    plot(ds.chnkrs, '-o', 'Color', cq)
    hold on
    quiver(ds.chnkrs, 'Color', cO)
    axis equal
    title("Capacitance Problem - Circles and Normals")
    % Plot the uks on surface
    figure()
    for i =1:n
        chnkri = ds.listChnkrs(i); 
        xOnSurface = reshape(chnkri.r, 2, chnkri.k*chnkri.nch);
        u_toUse = uk{i}(xOnSurface);
        scatter3(xOnSurface(1, :), xOnSurface(2, :), u_toUse, [], cO)
        hold on
    end
    title("Capacitance Problem - uk on surface")
    hold off
end


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
    nCol = nk(k + 1) - nk(k); % Number of columns for the submatrices
    chnkrk = ds.listChnkrs(k);
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nk(k) + 1;
    end_col = nk(k+1);
    for i=1:n
        nRow = nk(i + 1) - nk(i); % Number of rows for the submatrix
        start_row = nk(i) + 1;
        end_row = nk(i + 1);
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
        K(start_col:end_col, start_row:end_row) = submat;
    end
end



% Fill the RHS

rhs = zeros(Ntot, 1);

for i = 1:n
    start_row = nk(i) + 1;
    end_row = nk(i + 1);
    chnkri = ds.listChnkrs(i); 
    xOnSurface = reshape(chnkri.r, 2, chnkri.k*chnkri.nch);
    u_toUse = uk{i}(xOnSurface);
    rhs(start_row:end_row) = u_toUse;
end


% Solve for sigma, the unknown density

s = tic();
sigma = gmres(K, rhs, [], 1e-14, 100);
t2 = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with GMRES\n", t2);

if(plt)
    % Plot the sigma_ks on surface
    allXonSurface = reshape(ds.chnkrs.r, 2, ds.chnkrs.k*ds.chnkrs.nch);
    figure()
    scatter3(allXonSurface(1, :), allXonSurface(2, :), sigma, [], cq)
    title("Capacitance Problem - sigma on surface, GMRES")
end

% Compute the boundary integral of sigma over the discs to find q
q = zeros(n, 1);
opts.usesmooth = false;

for i=1:n
    chnkri = ds.listChnkrs(i);
    start_row = nk(i) + 1;
    end_row = nk(i + 1);
    sigma_i = sigma(start_row:end_row);
    q(i) = chunkerintegral(chnkri, sigma_i, opts);
    fprintf("%5.5e : charge at disk %1.0e with GMRES\n\n", q(i), i);
end


if(plt)
    % Plot the charges qk
    figure()
    for i=1:n
        [X,Y,Z] = ellipsoid( ctrs(1, i) , ctrs(2, i) , q(i) , ds.Rs(i) , ds.Rs(i) , 0 );
        s = surf(X,Y,Z,'FaceAlpha',0.5);
        s.EdgeColor = 'none';
        hold on
    end
    colorbar
    title("Capacitance Problem - Total charge on each disc, GMRES")
    hold off    

end

if( nargout > 2 )
    % We also need zztarg, xxtarg, yytarg
    % Find points off surface
    rmin = min(ds.chnkrs);
    rmin = rmin - outopt.fact*[1;1];
    rmax = max(ds.chnkrs);
    rmax = rmax + outopt.fact*[1;1];
    nplot = outopt.nplot;
    hx = (rmax(1)-rmin(1))/nplot;
    hy = (rmax(2)-rmin(2))/nplot;
    xtarg = linspace(rmin(1)+hx/2,rmax(1)-hx/2,nplot); 
    ytarg = linspace(rmin(2)+hy/2,rmax(2)-hy/2,nplot);
    [xxtarg,yytarg] = meshgrid(xtarg,ytarg);
    targets = zeros(2,length(xxtarg(:)));
    targets(1,:) = xxtarg(:); 
    targets(2,:) = yytarg(:);
    
    s = tic; 
    in = chunkerinterior(ds.chnkrs , targets ); 
    t3 = toc(s);
    fprintf('%5.2e s : time to find points off surface\n',t3)
    
    % Evaluate points off surface
    s = tic;
    Dsol_capacitance = chunkerkerneval(ds.chnkrs , DLplusSL,  sigma , targets(:,~in)); 
    t4 = toc(s);
    fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t4);
    
    % Add the off surface evaluations
    zztarg = nan(size(xxtarg));
    zztarg(~in) = Dsol_capacitance;
    
    if(plt)
        %%%%% Plot u off surface
        figure()
        h=surf(xxtarg,yytarg,zztarg);
        set(h,'EdgeColor','none')
        title("Capacitance Problem - Potential in the exterior")
        colorbar
    
    end
end



end