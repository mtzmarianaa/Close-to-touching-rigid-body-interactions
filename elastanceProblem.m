function [uk, sigma, nGMRES, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, solveType, plt, outopt)
% Given the description of the geometry of nCircles, solve the elastance
% problem
% IN: ds : a discs object with the given geometry of the discs
%        qk:  charges on the n given discs
%        solveType: 'full' solves the full system, 'precond' solves the
%                         precond full system, 'precondcomp' solves the precond compressed
%                         system
%        plt: boolean, wheather to plot or not
% OUT:  uk: u defined at the boundary
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


nDiscs = length(qk); % Number of discs
ctrs = ds.ctrs;

flagFunction = @(k, ds) dsc.flagnDisc(k, ds);
M = dsc.buildMelastance(ds, ds.listChnkrs, ds.nB, flagFunction);
kern = @(s,t) krns.DL_kern(s,t);
matOffSet = 0.5*eye(ds.chnkrs.npt) + M;

% Build according to preferences

if strcmp(solveType, 'full') || strcmp(solveType, 'precond')

    % Complete building the rhs 
    rhs = buildRHS_elastance(ds, ds.listChnkrs, ds.nB, flagFunction, qk);
else
    flagFunctionCoarse = @(k, ds) dsc.flagnDiscCoarse(k, ds);
    MCoarse = dsc.buildMelastance(ds, [ds.gamma0, ds.listCoarseGammas] , ds.nBCoarse, flagFunctionCoarse);
    matOffSetCoarse = 0.5*eye(ds.nBCoarse(end)) + MCoarse;

    % Build the rhs
    rhs = buildRHS_elastance(ds, [ds.gamma0, ds.listCoarseGammas] , ds.nBCoarse, flagFunctionCoarse, qk);
end


% See if we have the correct files for the interpolation
if strcmp(solveType, 'interprecondcomp')
    if isfile('../+prc/matInterpolant_Elastance.mat')
        load('../+prc/matInterpolant_Elastance.mat', 'matInterpolant_Elastance');
        matInterpolant = matInterpolant_Elastance;
    else
        % Meaning mat interpolant is not saved, but maybe we do have the
        % list of precomputed R
        if isfile('../+prc/listPrecomputedR_Elastance.mat')
            load('../+prc/listPrecomputedR_Elastance.mat', 'listPrecomputedR_Elastance');
            % With this we can build matInterpolant
            matInterpolant = rcip.buildInterp(listPrecomputedR_Elastance);
        else
            % We dont have mat interpolant and we dont have the list of
            % precomputed Rs
            geom0 = [];
            geom0.Rs = [0.75; 0.75];
            geom0.ctrs = [0  1.6; 0 0];
            pClose0 = [];
            pClose0(1).data = [0 2 1];
            pClose0(1).nClose = 1;
            pClose0(1).thetasReg = pi/6;
            pClose0(2).data = [pi, 1, 1];
            pClose0(2).nClose =1;
            pClose0(2).thetasReg = pi/6;
            pClose0(1).nBreakPoints = [10;10];
            pClose0(2).nBreakPoints = [10;10];
            if isfile('../+prc/listK22_invElastance.mat')
                load('../+prc/listK22_invElastance.mat', 'listK22_invElastance');
                [listPrecomputedR_Elastance, ~] = rcip.buildPrecomputedR_twoDiscs(geom0,  ...
                    pClose0, listK22_invElastance);
                matInterpolant = rcip.buildInterp(listPrecomputedR_Elastance);
            else
                listK22_invElastance = rcip.listK22_invElastance(geom0, pClose0);
                [listPrecomputedR_Elastance, ~] = rcip.buildPrecomputedR_twoDiscs(geom0,  ...
                    pClose0, listK22_invElastance);
                matInterpolant = rcip.buildInterp(listPrecomputedR_Elastance);
            end
        end
    end
end



% Solve according to preferences

if strcmp(solveType, 'full')
    [sigma, nGMRES] = dsc.solveFull(ds, rhs, kern, matOffSet);
elseif strcmp(solveType, 'precond')
    [sigma, nGMRES] = dsc.solveBlockPrecond(ds, rhs, kern, matOffSet);
elseif strcmp(solveType, 'precondcomp')
    [sigma, nGMRES] = dsc.solvePrecondComp(ds, rhs, kern, matOffSet, matOffSetCoarse);
elseif strcmp(solveType, 'interprecondcomp')
    [sigma, nGMRES] = dsc.solveInterpPrecond(ds, rhs, kern, matInterpolant, matOffSet, matOffSetCoarse);
end


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
    DL_kern = @(s,t) krns.DL_kern(s,t);
    SL_kern = @(s,t) krns.SL_kern(s,t);
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