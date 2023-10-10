function [uk, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, solveType, typeNodes, plt, outopt)
% *elastanceProblem* solves the elastance problem on non overlapping
% identical discs. Depending on the arguments it solves the problem with a
% different approach.
%
% Syntax: [uk, sigma, nGMRES] = elastanceProblem(ds, qk)
%              [uk, sigma, nGMRES] = elastanceProblem(ds, qk, solveType)
%              [uk, sigma, nGMRES] = elastanceProblem(ds, qk, solveType, typeNodes)
%              [uk, sigma, nGMRES] = elastanceProblem(ds, qk, solveType, typeNodes, plt)
%              [uk, sigma, nGMRES] = elastanceProblem(ds, qk, solveType, typeNodes, outopt)
%              [uk, sigma, nGMRES, tSolve] = elastanceProblem(ds, qk)
%              [uk, sigma, nGMRES, tSolve] = elastanceProblem(ds, qk, solveType)
%              [uk, sigma, nGMRES, tSolve] = elastanceProblem(ds, qk, solveType, typeNodes)
%              [uk, sigma, nGMRES, tSolve] = elastanceProblem(ds, qk, solveType, typeNodes, plt)
%              [uk, sigma, nGMRES, tSolve] = elastanceProblem(ds, qk, solveType, typeNodes, outopt)
%              [uk, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk)
%              [uk, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, solveType)
%              [uk, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, solveType, typeNodes)
%              [uk, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, solveType, typeNodes, plt)
%              [uk, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = elastanceProblem(ds, qk, solveType, typeNodes, outopt)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   qk - boundary information (charge on the discs)
%
% Optional input:
%   solveType - which type of solver to use 
%                                'full' solves the full system
%                                'precond' solves the preconditioned full system
%                                'precondcomp' solves the preconditioned compressed system
%                                'interprecondcomp'solves the
%                                preconditioned compressed system using
%                                interpolation
%   typeNodes - 'l' for Legendre nodes, 'logc' for log Chebyshev
%   plt - boolean if the method should render plots or not
%   outopt - output options for plotting
%
% Output:
%   uk - solution on the discs (organized by blocks)
%   sigma - solution density (organized by blocks)
%   nGMRES - number of GMRES iterations used to solve the system
%
% Optional output:
%   tSolve - time taken to assemble and solve the problem
%   zztarg - solution for plotting outside the discs
%   xxtarg - x coordinates for plotting outside the discs
%   yytarg - y coordinates for plotting outside the discs
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

if(nargin < 3)
    solveType = 'full';
end
solveType = lower(solveType);

if(nargin < 4)
    typeNodes = 'logc';
end
typeNodes = lower(typeNodes);

if( nargin < 5 )
    plt = false;
end

if( nargin < 6 )
    outopt = [];
    outopt.fact = 1.5;
    outopt.nplot = 750;
end

if ~isfield(outopt, 'fact')
    outopt.fact = 3;
end


if ~isfield(outopt, 'nplot')
    outopt.nplot = 250;
end

verbose = true;
if isfield(outopt, 'verbose')
    verbose = outopt.verbose;
end


nDiscs = length(qk); % Number of discs
ctrs = ds.ctrs;

flagFunction = @(k, ds) dsc.flagnDisc(k, ds);
M = dsc.elst.buildMelastance(ds, ds.listChnkrs, ds.nB, flagFunction);
kern = kernel('lap', 'd');
matOffSet = 0.5*eye(ds.chnkrs.npt) + M;

% Build according to preferences

if strcmp(solveType, 'full') || strcmp(solveType, 'precond')

    % Complete building the rhs 
    rhs = dsc.elst.buildRHS_elastance(ds, ds.listChnkrs, ds.nB, flagFunction, qk);
else
    flagFunctionCoarse = @(k, ds) dsc.flagnDiscCoarse(k, ds);
    MCoarse = dsc.elst.buildMelastance(ds, [ds.gamma0, ds.listCoarseGammas] , ds.nBCoarse, flagFunctionCoarse);
    matOffSetCoarse = 0.5*eye(ds.nBCoarse(end)) + MCoarse;

    % Build the rhs
    rhs = dsc.elst.buildRHS_elastance(ds, [ds.gamma0, ds.listCoarseGammas] , ds.nBCoarse, flagFunctionCoarse, qk);
    
end

% Save nu in the fine discretization if we need to plot
if(plt)
    [~, nu] = dsc.elst.buildRHS_elastance(ds, ds.listChnkrs, ds.nB, flagFunction, qk);
end

% See if we have the correct files for the interpolation
if strcmp(solveType, 'interprecondcomp')
    filePrecomputedR = strcat('../+prc/listPrecomputedR_Elastance', typeNodes ,'.mat');
    if isfile(filePrecomputedR)
        load(filePrecomputedR, 'listPrecomputedR_Elastance');
    else
        % We dont have mat interpolant and we dont have the list of
        % precomputed Rs
        fileListK22Inv = strcat('../+prc/listK22_invElastance.mat', typeNodes ,'.mat');
        geom0 = [];
        geom0.Rs = [0.75; 0.75];
        geom0.ctrs = [0  1.6; 0 0];
        geom0.nBreakPoints = [10; 10];
        pClose0 = [];
        pClose0(1).data = [0 2 1];
        pClose0(1).nClose = 1;
        pClose0(1).thetasReg = pi/6;
        pClose0(2).data = [pi, 1, 1];
        pClose0(2).nClose =1;
        pClose0(2).thetasReg = pi/6;
        pClose0(1).nBreakPoints = [10;10];
        pClose0(2).nBreakPoints = [10;10];
        if isfile(fileListK22Inv)
            load(fileListK22Inv, 'listK22_invElastance');
            listPrecomputedR_Elastance = rcip.buildPrecomputedR_twoDiscs(geom0,  ...
                pClose0, listK22_invElastance, typeNodes);
            save(filePrecomputedR, 'listPrecomputedR_Elastance');
        else
            listK22_invElastance = dsc.elst.listK22_invElastance(geom0, pClose0, typeNodes);
            save(fileListK22Inv, 'listK22_invElastance');
            listPrecomputedR_Elastance = rcip.buildPrecomputedR_twoDiscs(geom0,  ...
                pClose0, listK22_invElastance, typeNodes);
            save(filePrecomputedR, 'listPrecomputedR_Elastance');
        end
    end
end

% Solve according to preferences

if strcmp(solveType, 'full')
    [sigma, nGMRES, tS] = dsc.solveFull(ds, rhs, kern, matOffSet, verbose);
elseif strcmp(solveType, 'precond')
    [sigma, nGMRES, tS] = dsc.solveBlockPrecond(ds, rhs, kern, matOffSet, verbose);
elseif strcmp(solveType, 'precondcomp')
    [sigma, nGMRES, tS] = dsc.solvePrecondComp(ds, rhs, kern, matOffSet, matOffSetCoarse, 0, 0, 0, verbose);
elseif strcmp(solveType, 'interprecondcomp')
    [sigma, nGMRES, tS] = dsc.solveInterpPrecond(ds, rhs, kern, listPrecomputedR_Elastance, matOffSet, matOffSetCoarse, 0, verbose);
end


% Compute u for on surface points

uk = -M*sigma;


if(nargout > 3)
    tSolve = tS;
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


if( nargout > 4 )
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
    if(verbose)
    fprintf('%5.2e s : time to find points off surface\n',t3)
    end
    
    % Evaluate points off surface
    DL_kern = @(s,t) krns.DL_kern(s,t);
    SL_kern = @(s,t) krns.SL_kern(s,t);
    s = tic;
    Dsol_elastance = chunkerkerneval(ds.chnkrs , DL_kern,  sigma , targets(:,~in));
    Dsol_elastance = Dsol_elastance + chunkerkerneval(ds.chnkrs, SL_kern, nu, targets(:, ~in));
    t4 = toc(s);
    if(verbose)
    fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t4);
    end
    
    % Add the off surface evaluations
    zztarg = nan(size(xxtarg));
    zztarg(~in) = Dsol_elastance;
    
    if(plt)
        %%%%% Plot u off surface
        figure()
        h=surf(xxtarg,yytarg,zztarg);
        set(h,'EdgeColor','none')
        title("Elastance problem - Potential in the exterior")
        colorbar
        view(2)
        axis equal
        grid off
        ylim([-1.5 1.5])
    
    end
end




end