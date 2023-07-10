%%% Capacitance problem test solution density vs angle

addpaths_loc();
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
nBreakPoints = [20; 20];

%%%%%%%%%% UNIFORM PANELS
% Define points on surface
geom = [];
geom.Rs = Rs;
geom.nBreakPoints = nBreakPoints;

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
points1 = reshape( ds2.listChnkrs(1).r, 2, [] );
points2 = reshape( ds2.listChnkrs(2).r, 2, [] );
thetas1 = acos(  (points1(1, :) - ctrs(1, 1) )./( Rs(1)  ) );
thetas2 = acos(  (points2(1, :) - ctrs(1, 2) )./( Rs(2) ) );


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





