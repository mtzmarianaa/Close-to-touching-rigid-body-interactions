function [R, G] = buildR(chunkerCoarse, chunkerFine, K22_inv, P,  kernel)
% *buildR* build the matrix for RCIP 
% $$R = W_c^{-1} P' W_f K_{22}^{-1} P. $$
%
% Syntax: R = buildR(chunkerCoarse, chunkerFine)
%              R = buildR(chunkerCoarse, chunkerFine, K22_inv, P,  kernel)
%              [R, G] = buildR(chunkerCoarse, chunkerFine)
%              [R, G] = buildR(chunkerCoarse, chunkerFine, K22_inv, P,  kernel)
%
% Input:
%   chunkerCoarse - chunker object describing the coarse discretization
%   chunkerFine - chunker object describing the fine discretization
%
% Optional input:
%   K22_inv - inverse for the close-to-touching region
%   P - prolongation matrix
%   kernel - kernel to use in the BIE
%
% Output:
%   R - $R = W_c^{-1} P' W_f K_{22}^{-1} P. $
%   G - $G = W_c^{-1} P' W_f. $
%
% author: Mariana Martinez (mariana.martinez.aguilar@gmail.com)

% We need Wc, Wf, P, K22
Wc = weights(chunkerCoarse);
Wc = 1./(Wc(:));
Wc = spdiags(Wc);
nc = chunkerCoarse.npt;
Wc = speye(nc).*Wc;

Wf = weights(chunkerFine);
Wf = Wf(:);
Wf = spdiags(Wf);
nf = chunkerFine.npt;
Wf = speye(nf).*Wf;

nRef = floor( chunkerFine.nch/4 - 2  );
nRef = max(0, nRef);

if(nargin < 4)
    P = rcip.prol_dyadic(chunkerCoarse.k, nRef);
    P = blkdiag(P, P);
end

if(nargin < 3)
    K22 = chunkermat(chunkerFine, kernel);
    K22_inv = inv(K22);
end

% Build R and G if needed
 if( nargout < 2) 
     R = Wc*P'*Wf*K22_inv*P;
 else
     G = Wc*P'*Wf;
     R = G*K22_inv*P;
 end


end