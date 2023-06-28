%%% Capacitance problem test solution density vs angle


clear all
close all
clc

u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);

nTest = 5;

uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 1.5005 ;0 0]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [10; 10];


% Define points on surface
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
geom.saveAngles = true;

k = 16;
x = lege.exps(k);
sigma1 = zeros(nTest, k*(nBreakPoints(1)-1));
sigma2 = zeros(nTest, k*(nBreakPoints(2) -1) );

xcoordCtr2 = linspace(1.5005, 5, nTest ); % Variation in the coordinate ctrs(1,2), moving the discs closer

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    geom.ctrs = ctrs;
    ds = discs(geom);
    [q, sigma] = capacitanceProblem(ds, uk);
    % Fill the matricies with computed densities
    sigma1(i, :) = sigma(1:k*(nBreakPoints(1)-1) ) ; 
    sigma2(i, :) = sigma(( k*(nBreakPoints(1)-1) + 1 ):end);
end

% Put everything in order
k = ds.listChnkrs.k;
x = lege.exps(k); % Get the legendre nodes

% Get the angles
thetas1 = zeros(1, k*nBreakPoints(1));
thetas2 = zeros(1, k*nBreakPoints(2));

% Add the breakpoints
thetas1(1:k:end) = ds.thetas(1).breakpoints;
thetas2(1:k:end) = ds.thetas(2).breakpoints;

% Add the inbetween nodes


% Plot density vs angle for both circles


