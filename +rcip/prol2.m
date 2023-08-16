function submat = prol2( x )
% *prol2* For x, n Legendre nodes, builds the prolongation matrix that
% interpolates data from this discretization to a finer one with 2*n
% Gauss-Legendre nodes. Done by solving two Vandermonde systems.
%
% Syntax: submat = prol2( x )
%
% Input:
%   x - Gauss-Legendre nodes
%
% Output:
%   submat - prolongation matrix
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

n = length(x); % Get number of Gauss-Legendre nodes on the coarse discretization
% Build the Vandermonde matrices
x_fine = [x-1; x+1]; % Legendre nodes in the finer discretization
x_fine = x_fine./2;
C = rcip.shortVandermonde(x, n);
F = rcip.shortVandermonde(x_fine, n);
submat = F/C;
end

