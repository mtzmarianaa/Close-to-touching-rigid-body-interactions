% Capacitance problem (2 discs first)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
format long
% For a point charge
x0 = 0.75/(2*sqrt(2))*[1; 1];

% u1 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
% u2 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);


u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);

uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 1.5005 ;0 0]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];


% Define points on surface
geom = [];
geom.ctrs = ctrs;
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
ds = discs(geom);


[q, sigma, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, true);



% % Plot errors off surface
% targets = zeros(2,length(xxtarg(:)));
% targets(1,:) = xxtarg(:); 
% targets(2,:) = yytarg(:);
% trueTarg = nan(size(zztarg));
% in = chunkerinterior(ds.chnkrs , targets ); 
% trueTarg(~in) = u1(targets(:,~in));
% relErr = abs(zztarg - trueTarg)./abs(trueTarg);
% 
% figure()
% h=surf(xxtarg, yytarg, relErr );
% set(h,'EdgeColor','none')
% title("Capacitance Problem - Relative Error in Potential in the exterior")
% colorbar
% 
% 
% figure()
% x = [min(xxtarg(:)), max(xxtarg(:))];
% y = [min(yytarg(:)), max(yytarg(:))];
% imagesc(x, y, relErr)
% title("Capacitance Problem - Relative Error in Potential in the exterior")
% colorbar

