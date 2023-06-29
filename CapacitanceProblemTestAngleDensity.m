%%% Capacitance problem test solution density vs angle


clear all
close all
clc

u1 = @(x) 0*x(1, :);
u2 = @(x) 1+0*x(1, :);

nTest = 5;

uk = {u1, u2}; % Functions uk, u on the boundary of the k-th circle
ctrs = [0 1.500005 ;0 0]; % Centers of the circles
Rs = [0.75; 0.75]; % Radi of the circles
n = length(uk);
nBreakPoints = [12; 12];

%%%%%%%%%% UNIFORM PANELS
% Define points on surface
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
geom.saveAngles = true;

k = 16;
x = lege.exps(k);
sigma1 = zeros(nTest, k*(nBreakPoints(1)-1));
sigma2 = zeros(nTest, k*(nBreakPoints(2) -1) );

xcoordCtr2 = linspace(1.500005, 1.75, nTest ); % Variation in the coordinate ctrs(1,2), moving the discs closer

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    geom.ctrs = ctrs;
    ds2 = discs(geom);
    [q, sigma] = capacitanceProblem(ds2, uk);
    % Fill the matricies with computed densities
    sigma1(i, :) = sigma(1:k*(nBreakPoints(1)-1) ) ; 
    sigma2(i, :) = sigma(( k*(nBreakPoints(1)-1) + 1 ):end);
end

% Put everything in order
k = ds2.listChnkrs.k;
x = lege.exps(k); % Get the legendre nodes


% Get the angles
thetas1 = chunksInParameter(ds2.listChnkrs(1), ds2.thetas(1).breakpoints);
thetas2 = chunksInParameter(ds2.listChnkrs(2), ds2.thetas(2).breakpoints);


% Plot density vs angle for both circles
% For circle1
figure()
leg = string(xcoordCtr2 - 1.5);
plot(thetas1(:), sigma1, '-o')
title("Solution density and angle, disc 1")
xlabel("Angle")
ylabel("\sigma")
legend(leg)


% For circle2
figure()
leg = string(xcoordCtr2 - 1.5);
plot(thetas2(:), sigma2, '-o')
title("Solution density and angle, disc 2")
xlabel("Angle")
ylabel("\sigma")
legend(leg)




%%%%%%%%%% NON UNIFORM PANELS
% Define points on surface
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;
geom.saveAngles = true;
geom.ctrs = ctrs;

pClose = [];
pClose(1).thetas = [0.0];
pClose(1).nClose = 1;
pClose(1).discClose = 2;
pClose(1).pRef = [pi];
pClose(2).thetas = [pi];
pClose(2).nClose = 1;
pClose(2).discClose = 1;
pClose(2).pRef = [0.0];

ds2 = discs(geom, pClose);

npt1 = ds2.listChnkrs(1).npt;
npt2 = ds2.listChnkrs(2).npt;

k = 16;
x = lege.exps(k);
sigma1 = zeros(nTest, npt1 );
sigma2 = zeros(nTest, npt2 );

xcoordCtr2 = linspace(1.500005, 1.75, nTest ); % Variation in the coordinate ctrs(1,2), moving the discs closer

for i=1:nTest
    ctrs(1,2) = xcoordCtr2(i);
    geom.ctrs = ctrs;
    ds2 = discs(geom, pClose);
    [q, sigma] = capacitanceProblem(ds2, uk);
    % Fill the matricies with computed densities
    sigma1(i, :) = sigma(1:npt1 ) ; 
    sigma2(i, :) = sigma(( npt1 + 1 ):end);
end

% Put everything in order
k = ds2.listChnkrs.k;
x = lege.exps(k); % Get the legendre nodes


% Get the angles
thetas1 = chunksInParameter(ds2.listChnkrs(1), ds2.thetas(1).breakpoints);
thetas2 = chunksInParameter(ds2.listChnkrs(2), ds2.thetas(2).breakpoints);


% Plot density vs angle for both circles
% For circle1
figure()
leg = string(xcoordCtr2 - 1.5);
plot(thetas1(:), sigma1, '-o')
title("Solution density and angle, disc 1")
xlabel("Angle")
ylabel("\sigma")
legend(leg)


% For circle2
figure()
leg = string(xcoordCtr2 - 1.5);
plot(thetas2(:), sigma2, '-o')
title("Solution density and angle, disc 2")
xlabel("Angle")
ylabel("\sigma")
legend(leg)

