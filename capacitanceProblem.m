function [q, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, solveType, typeNodes, plt, outopt)
% *capacitanceProblem* solves the capacitance problem on non overlapping
% identical discs. Depending on the arguments it solves the problem with a
% different approach.
%
% Syntax: [q, sigma, nGMRES] = capacitanceProblem(ds, uk)
%              [q, sigma, nGMRES] = capacitanceProblem(ds, uk, solveType)
%              [q, sigma, nGMRES] = capacitanceProblem(ds, uk, solveType, plt)
%              [q, sigma, nGMRES] = capacitanceProblem(ds, uk, solveType, outopt)
%              [q, sigma, nGMRES, tSolve] = capacitanceProblem(ds, uk)
%              [q, sigma, nGMRES, tSolve] = capacitanceProblem(ds, uk, solveType)
%              [q, sigma, nGMRES, tSolve] = capacitanceProblem(ds, uk, solveType, typeNodes)
%              [q, sigma, nGMRES, tSolve] = capacitanceProblem(ds, uk, solveType, typeNodes, plt)
%              [q, sigma, nGMRES, tSolve] = capacitanceProblem(ds, uk, solveType, typeNodes, outopt)
%              [q, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk)
%              [q, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, solveType)
%              [q, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, solveType, typeNodes)
%              [q, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, solveType, typeNodes, plt)
%              [q, sigma, nGMRES, tSolve, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, solveType, typeNodes, outopt)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   uk - boundary information
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
%   q - charges on each of the discs (constant per disc)
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

n = length(uk);
ctrs = ds.ctrs;

kern = kernel('lap', 'c', [1.0, 1.0]);
matOffSet = 0.5*eye(ds.chnkrs.npt);

% Build according to preferences
if strcmp(solveType, 'full') || strcmp(solveType, 'precond')
    % Build the rhs
    flagFunction = @(k, ds) dsc.flagnDisc(k, ds);
    rhs = dsc.capc.buildRHS_capacitance(ds, ds.listChnkrs, ds.nB, flagFunction, uk);
    
else
    % Build the rhs
    flagFunction = @(k, ds) dsc.flagnDiscCoarse(k, ds);
    rhs = dsc.capc.buildRHS_capacitance(ds, [ds.gamma0, ds.listCoarseGammas], ds.nBCoarse, flagFunction, uk);
    
    % Solve the system
    matOffSetCoarse = 0.5*eye(ds.nBCoarse(end));
end

% See if we have the correct files for the interpolation
if strcmp(solveType, 'interprecondcomp')
    filePrecomputedR = strcat('../+prc/listPrecomputedR_Capacitance', typeNodes ,'.mat');
    if isfile(filePrecomputedR)
        load(filePrecomputedR, 'listPrecomputedR_Capacitance');
    else
        % We dont have mat interpolant and we dont have the list of
        % precomputed Rs
        fileListK22Inv = strcat('../+prc/listK22_invCapacitance.mat', typeNodes ,'.mat');
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
            load(fileListK22Inv, 'listK22_invCapacitance');
            listPrecomputedR_Capacitance = rcip.buildPrecomputedR_twoDiscs(geom0,  ...
                pClose0, listK22_invCapacitance, typeNodes);
            save(filePrecomputedR, 'listPrecomputedR_Capacitance');
        else
            listK22_invCapacitance = dsc.capc.listK22_invCapacitance(geom0, pClose0, typeNodes);
            save(fileListK22Inv, 'listK22_invCapacitance');
            listPrecomputedR_Capacitance = rcip.buildPrecomputedR_twoDiscs(geom0,  ...
                pClose0, listK22_invCapacitance, typeNodes);
            save(filePrecomputedR, 'listPrecomputedR_Capacitance');
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
    [sigma, nGMRES, tS] = dsc.solveInterpPrecond(ds, rhs, kern, listPrecomputedR_Capacitance, typeNodes, matOffSet, matOffSetCoarse, 0, verbose);
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
    if(verbose)
    fprintf("%5.5e : charge at disk %1.0e with GMRES\n\n", q(i), i);
    end
end

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
    title("Capacitance Problem - sigma on surface")


    % Plot the charges qk
    maxQ = round( max(abs(q)) );
    figure()
    for i=1:n
        [X,Y,Z] = ellipsoid( ctrs(1, i) , ctrs(2, i) , q(i) , ds.Rs(i) , ds.Rs(i) , 0 , 50);
        s = surf(X,Y,Z,'FaceAlpha',0.9);
        s.EdgeColor = 'none';
        annotation = ['q_{', num2str(i), '}=', num2str(q(i))];
        text(ctrs(1, i) , ctrs(2, i) , q(i), annotation)
        hold on
    end
    colorbar
    view(2)
    axis equal
    clim([-maxQ, maxQ]);
    grid off
    title("Capacitance Problem - Total charge on each disc")
    ylim([-1.5 1.5])
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
    DLplusSL = kernel('lap', 'c', [1.0, 1.0]);
    s = tic;
    %Dsol_capacitance = DLplusSL.fmm(1e-12, ds.chnkrs, targets(:,~in), sigma);
    Dsol_capacitance = chunkerkerneval(ds.chnkrs , DLplusSL,  sigma , targets(:,~in)); 
    t4 = toc(s);
    if(verbose)
    fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t4);
    end
    
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