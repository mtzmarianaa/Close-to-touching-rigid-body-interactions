function [q, sigma, nGMRES, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, plt, outopt)
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
nChunkers = length(ds.listChnkrs);
ctrs = ds.ctrs;
Rs = ds.Rs;
nBreakPoints = ds.nBreakPoints;


Ntot = ds.chnkrs.npt;
nB = ds.nB;


% BUILD THE KERNEL FUNCTIONS
SL_kern = @(s,t) chnk.lap2d.kern(s, t, 's');
DL_kern = @(s,t) chnk.lap2d.kern(s, t, 'd');
DLplusSL = @(s,t) DL_kern(s,t) + SL_kern(s,t);
ID_kern = @(s,t)  speye(size(t.r, 2)*size(t.r, 3), size(s.r, 2) *size(s.r, 3));
DLplusSLplusHI = @(s,t) 0.5*ID_kern(s,t) + DLplusSL;


%%%%%%%%%%%
%%%%%% OPTION A: BUILD THE MATRIX DIRECTLY WITH CHUNKERMAT
% Build the matrix
% K = chunkermat( ds.chnkrs, DLplusSL);
% K = K + 0.5*eye(size(K));


%%%%%%%%% OPTION B: BUILD THE MATRIX BLOCK BY BLOCK
% Initialize the matrix
K = zeros(Ntot);
opts = [];
opts2 = [];
opts2.adaptive_correction = true;


% Proceed to fill in the matrix 
for k=1:nChunkers
    % fill column by column
    nCol = nB(k + 1) - nB(k); % Number of columns for the submatrices
    chnkrk = ds.listChnkrs(k);
    targ = reshape(chnkrk.r, 2 , chnkrk.k*chnkrk.nch);
    start_col = nB(k) + 1;
    end_col = nB(k+1);
    for i=1:nChunkers
        nRow = nB(i + 1) - nB(i); % Number of rows for the submatrix
        start_row = nB(i) + 1;
        end_row = nB(i + 1);
        chnkri = ds.listChnkrs(i); % chunker we are working with
        % See if we have to do an off boundary or on boundary eval
        if(i == k)
            % on boundary
            submat = chunkermat(chnkri, DLplusSL, opts2) + 0.5*eye(nRow);
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
    % Use the flag
    flag = logical( dsc.flagnDisc(i, ds) );
    % Convert that to points
    flag_points = repmat(flag, 1, ds.chnkrs.k);
    flag_points = flag_points';
    flag_points = logical( flag_points(:) );
    % Get the points
    xOnSurface = reshape(ds.chnkrs.r(:, :, flag ), 2, []);
    u_toUse = uk{i}(xOnSurface);
    rhs(flag_points) = u_toUse';
end


% Solve for sigma, the unknown density
s = tic();
if( nargout > 2 )
    [sigma, ~, ~, nGMRES] = gmres(K, rhs, [], 1e-14);
else
    sigma = gmres(K, rhs, [], 1e-14, 1500);
end
t2 = toc(s);
fprintf("%5.2e s :time taken to solve the linear system with GMRES\n", t2);


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