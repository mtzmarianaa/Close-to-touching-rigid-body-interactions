% Test the capacitance problem without adaptivity and by moving the discs
% closer together

clear all
close all
clc

nTest = 100;
plt = false;
% For a point charge
x0 = 0.75/(2*sqrt(2))*[1; 1];

u1 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);
u2 = @(x) log(1./vecnorm( bsxfun(@minus, x, x0) ))/(2.0*pi);


uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 5 ;0 0]; % Centers of the circles, just going to change ctrs(1, 2)
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];


xcoordCtr2 = linspace(1.5005, 5, nTest ); % Variation in the coordinate ctrs(1,2), moving the discs closer
absErrs = zeros(nTest, 1);
relErrs = zeros(nTest, 1);
% Iterate

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    % Define points on surface
    geom = [];
    geom.ctrs = ctrs;
    geom.Rs = Rs;
    geom.nBreakPoints = nBreakPoints;
    ds = discs(geom);
    % Solve
    [q, sigma, zztarg, xxtarg, yytarg] = capacitanceProblem(ds, uk, false);
    % Compute the absolute and relative errors off surface
    targets = zeros(2,length(xxtarg(:)));
    targets(1,:) = xxtarg(:); 
    targets(2,:) = yytarg(:);
    trueTarg = nan(size(zztarg));
    in = chunkerinterior(ds.chnkrs , targets ); 
    trueTarg(~in) = u1(targets(:,~in) );
    relErr = abs(zztarg - trueTarg)./abs(trueTarg);
    absErrs(i) = norm(abs(zztarg(~in) - trueTarg(~in)));
    relErrs(i) = norm(relErr(~in));

    % Plot if necessary
    if(plt)
%     figure()
%     h=surf(xxtarg, yytarg, relErr );
%     set(h,'EdgeColor','none')
%     title("Capacitance Problem - Relative Error in Potential in the exterior")
%     colorbar
    
    figure()
    x = [min(xxtarg(:)), max(xxtarg(:))];
    y = [min(yytarg(:)), max(yytarg(:))];
    imagesc(x, y, relErr)
    title("Capacitance Problem - Relative Error in Potential in the exterior")
    colorbar
    end
end


% Plot the errors vs distance of circles

distCirc = xcoordCtr2 - 1.5;
figure()
plot( distCirc, absErrs, '-o' )
title("Distance between circles and absolute error potential in exterior")
xlabel("Distance between circles")
ylabel("Absolute error")

figure()
plot( distCirc, relErrs, '-o' )
title("Distance between circles and relative error potential in exterior")
xlabel("Distance between circles")
ylabel("Relative error")

