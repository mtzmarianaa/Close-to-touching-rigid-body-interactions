function [uk, sigma, nGMRES, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, plt, outopt)
% Given the description of the geometry of nCircles, solve the elastance
% problem
% IN: ds : a discs object with the given geometry of the discs
%        qk:  charges on the n given discs
%        plot: boolean, wheather to plot or not
% OUT:  uk: u defined at the boundary
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


nDiscs = length(qk); % Number of discs
nChunkers = length(ds.listChnkrs); % Number of chunker objects
ctrs = ds.ctrs;
Rs = ds.Rs;
nBreakPoints = ds.nBreakPoints;

Ntot = ds.chnkrs.npt; % Number of points in the whole problem
nB = ds.nB;


% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');
DLplusSL = @(s,t) DL_kern(s,t) + SL_kern(s,t);
ID_kern = @(s,t)  speye(size(t.r, 2)*size(t.r, 3), size(s.r, 2) *size(s.r, 3));


% Define the points on surface
xOnSurface = reshape(ds.chnkrs.r, 2, ds.chnkrs.k*ds.chnkrs.nch);


% Build nu (this is used for the RHS)
nu = ones(ds.chnkrs.npt, 1);
for i= 1:nDiscs
    % Use the flag
    flag = logical( dsc.flagnDisc(i, ds) );
    % Convert that to points
    flag_points = repmat(flag, 1, ds.chnkrs.k);
    flag_points = flag_points';
    flag_points = logical( flag_points(:) );
    f = ones(ds.chnkrs.npt, 1);
    f( ~flag_points) = 0; % Cancel out contribution from other discs
    % Get the perimeter of omegak
    per = chunkerintegral(ds.chnkrs, f, []);
    nu(flag_points) = qk(i)/per;
end

% Build the RHS matrix, the I/2 + DL matrix and the M matrix
IDplusDLMat = zeros(Ntot);
rhs_chunkermat = zeros(Ntot);
opts = [];

for k=1:nChunkers
    % Fill column by column
    nCol = nB(k+1) - nB(k); % Number of columns for the submatrices
    chnkrk = ds.listChnkrs(k); % Current chunker (i.e. current gamma)
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nB(k) + 1;
    end_col = nB(k+1);
    for i=1:nChunkers
        nRow = nB(i + 1) - nB(i); % Number of rows for the submatrix
        start_row = nB(i) + 1;
        end_row = nB(i+1) ;
        chnkri = ds.listChnkrs(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i == k)
            % on boundary
            submat = chunkermat(chnkri, DL_kern) ;
            submat = submat + 0.5*eye( size(submat) );
            submat_rhs = chunkermat(chnkri, SL_kern);
        else
            % off boundary
            submat = chunkerkernevalmat(chnkri, DL_kern, targ, opts);
            submat_rhs = chunkerkernevalmat(chnkri, SL_kern, targ, opts);
        end
        % Add this submatrix to K
        IDplusDLMat(start_col:end_col, start_row:end_row) = submat;
        rhs_chunkermat(start_col:end_col, start_row:end_row) = submat_rhs;
    end
end


% Build the block diagonal matrix M (note that here we don't care about the
% order of the discs, just the order of the chunker objects - order of
% either discs (if no close information) or the order of gammas (if close
% information)
M = zeros(Ntot);
% We need to "filter" out the weights from the whole chnkrs according to
% which disc they are found in
all_weights = weights( ds.chnkrs );
all_weights = reshape(all_weights, ds.chnkrs.npt, [])';

for i=1:nChunkers
    start_index = nB(i) + 1;
    end_index = nB(i+1);
    % Now filter according to discs
    for k=1:nDiscs
        flag = logical( dsc.flagnDisc(k, ds) );
        % Convert that to points
        flag_points = repmat(flag, 1, ds.chnkrs.k);
        flag_points = flag_points';
        flag_points = logical( flag_points(:) );
        flag_points_here = false( 1, Ntot);
        flag_points_here(start_index:end_index) = flag_points(start_index:end_index);
        % Fill the matrix
        weights_here = zeros( 1, Ntot );
        weights_here(flag_points') = all_weights(flag_points');
        weights_here = repmat(weights_here, sum(flag_points(start_index:end_index)), 1);
        M( flag_points_here, :) = weights_here;
    end
end



% Complete building the rhs 
%rhs_chunkermat = chunkermat(ds.chnkrs, SL_kern,  []);
rhs = -rhs_chunkermat*nu;


% Build the system
% IDplusDLMat = chunkermat(ds.chnkrs, DL_kern, opts);
% IDplusDLMat = IDplusDLMat + 0.5*speye(size(IDplusDLMat));
K = IDplusDLMat + M;

% Solve for sigma, unknown density
s = tic();
if( nargout > 2 )
    [sigma, ~, ~, nGMRES] = gmres(K, rhs, [], 1e-14, 100);
else
    sigma = gmres(K, rhs, [], 1e-14, 100);
end
t = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with Matlab's backslash\n", t);


% Compute u for on surface points

uk = -M*sigma;


% Plot if necessary

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
    title("Elastance Problem - Circles and Normals")
    hold off

    % Plot the carges
    figure()
    for i=1:nDiscs
        [X,Y,Z] = ellipsoid( ctrs(1, i) , ctrs(2, i) , qk(i) , ds.Rs(i) , ds.Rs(i) , 0 );
        s = surf(X,Y,Z,'FaceAlpha',0.5);
        s.EdgeColor = 'none';
        hold on
    end
    colorbar
    title("Elastance Problem - Total charge on each disc")
    hold off

    % Plot the sigma_ks on surface
    allXonSurface = reshape(ds.chnkrs.r, 2, ds.chnkrs.k*ds.chnkrs.nch);
    figure()
    scatter3(allXonSurface(1, :), allXonSurface(2, :), sigma, [], cq)
    title("Elastance Problem - sigma on surface, GMRES")

    % Plot the uks on surface - we separate everything to plot with
    % different colors
    figure()
    for i =1:nDiscs
        % Use the flag
        flag = logical( dsc.flagnDisc(i, ds) );
        % Convert that to points
        flag_points = repmat(flag, 1, ds.chnkrs.k);
        flag_points = flag_points';
        flag_points = logical( flag_points(:) );
        % Get the points
        xOnSurface = reshape(ds.chnkrs.r(:, :, flag ), 2, []);
        u_toUse = uk(flag_points);
        scatter3(xOnSurface(1, :), xOnSurface(2, :), u_toUse, [], cO)
        hold on
    end
    title("Elastance Problem - uk on surface, GMRES")
    hold off

end


if( nargout > 3 )
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
    Dsol_elastance = chunkerkerneval(ds.chnkrs , DL_kern,  sigma , targets(:,~in));
    Dsol_elastance = Dsol_elastance + chunkerkerneval(ds.chnkrs, SL_kern, nu, targets(:, ~in));
    t4 = toc(s);
    fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t4);
    
    % Add the off surface evaluations
    zztarg = nan(size(xxtarg));
    zztarg(~in) = Dsol_elastance;
    
    if(plt)
        %%%%% Plot u off surface
        figure()
        h=surf(xxtarg,yytarg,zztarg);
        set(h,'EdgeColor','none')
        title("Elastance Problem - Potential in the exterior")
        colorbar
    
    end
end




end