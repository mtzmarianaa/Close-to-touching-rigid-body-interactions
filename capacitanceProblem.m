function [q, sigma, nGMRES, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, solveType, plt, outopt)
% Given the description of the geometry of nCircles solve the capacitance
% problem. 
% IN: ds : a discs object with the given geometry of the discs
%        uk:  u defined at the boundary
%        solveType: 'full' solves the full system, 'precond' solves the
%                         precond full system, 'precondcomp' solves the precond compressed
%                         system
%        plt: boolean, wheather to plot or not
% OUT:  q: charges on the n given discs
%           sigma: on bundary density

if(nargin < 3)
    solveType = 'full';
end
solveType = lower(solveType);

if( nargin < 4 )
    plt = false;
end
if( nargin < 5 )
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

kern = @(s,t) krns.DLplusSL(s,t);
matOffSet = 0.5*eye(ds.chnkrs.npt);

% Build according to preferences
if strcmp(solveType, 'full') || strcmp(solveType, 'precond')
    % Build the rhs
    flagFunction = @(k, ds) dsc.flagnDisc(k, ds);
    rhs = buildRHS_capacitance(ds, ds.listChnkrs, ds.nB, flagFunction, uk);
    
else
    % Build the rhs
    flagFunction = @(k, ds) dsc.flagnDiscCoarse(k, ds);
    rhs = buildRHS_capacitance(ds, [ds.gamma0, ds.listCoarseGammas], ds.nBCoarse, flagFunction, uk);
    
    % Solve the system
    matOffSetCoarse = 0.5*eye(ds.nBCoarse(end));
end

% Solve according to preferences
if strcmp(solveType, 'full')
    [sigma, nGMRES] = dsc.solveFull(ds, rhs, kern, matOffSet);
elseif strcmp(solveType, 'precond')
    [sigma, nGMRES] = dsc.solveBlockPrecond(ds, rhs, kern, matOffSet);
else
    [sigma, nGMRES] = dsc.solvePrecondComp(ds, rhs, kern, matOffSet, matOffSetCoarse);
end


% Compute the boundary integral of sigma over the discs to find q
q = zeros(n, 1);
opts.usesmooth = false;

for i=1:n
    % Use the flag
    flag = logical( dsc.flagnDisc(i, ds) );
    % Convert that to points
    flag_points = repmat(flag, 1, ds.chnkrs.k);
    flag_points = flag_points';
    flag_points = logical( flag_points(:) );
    % Get the correct sigma_i
    sigma_i = sigma;
    sigma_i(~flag_points) = 0; % Set all to zero if not in disc, so we just get values from the discs
    q(i) = chunkerintegral(ds.chnkrs, sigma_i, opts);
    fprintf("%5.5e : charge at disk %1.0e with GMRES\n\n", q(i), i);
end



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
    title("Capacitance Problem - Circles and Normals")

    % Plot the uks on surface
    figure()
    for i =1:n
        % Use the flag
        flag = logical( dsc.flagnDisc(i, ds) );
        % Get the points
        xOnSurface = reshape(ds.chnkrs.r(:, :, flag ), 2, []);
        u_toUse = uk{i}(xOnSurface);
        scatter3(xOnSurface(1, :), xOnSurface(2, :), u_toUse, [], cO)
        hold on
    end
    title("Capacitance Problem - uk on surface")
    hold off

    % Plot the sigma_ks on surface
    allXonSurface = reshape(ds.chnkrs.r, 2, ds.chnkrs.k*ds.chnkrs.nch);
    figure()
    scatter3(allXonSurface(1, :), allXonSurface(2, :), sigma, [], cq)
    title("Capacitance Problem - sigma on surface, GMRES")


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