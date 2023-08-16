function listKs_inv = computeKii_inv(ds, kernel, matOffSet)
% *computeKii_inv* computes the list of inverses of the
% close-to-touching regions (useful for RCIP).
%
% Syntax: listKs_inv = computeKii_inv(ds, kernel, matOffSet)
%
% Input:
%   ds - discs object, has all the geometric properties of the collection
%          of non overlapping discs, their close-to-touching regions and their far
%          regions.
%   kernel - kernel object (from chunkie) or function handle definind the
%                kernel to use
%   matOffSet - matrix defining the integral operator which is not a kernel
%                        (has to be a matrix)
%
% Output:
%   listKs_inv - list of inverses for the close-to-touching region,
%                      $K_{22}^{-1}, K_{33}^{-1}, ..., K_{NN}^{-1}$ (Stores
%                      them in a cell array)
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

% Determine if matrix is given, if not initialize it with zeros
if( nargin < 3)
    matOffSet = zeros(ds.chnkrs.npt);
end

opts2 = [];
opts2.adaptive_correction = true;

% Compute the inverses
listKs_inv = cell(1, ds.nGammas);
for i=1:ds.nGammas
    Kii = chunkermat(ds.listGammas(i), kernel, opts2);
    Kii = Kii + matOffSet( (ds.nB(i+1)+1):ds.nB(i+2), (ds.nB(i+1)+1):ds.nB(i+2) );
    listKs_inv{i} = inv(Kii);
end


end