function listKs_inv = computeKii_inv(ds, kernel, matOffSet)
% Builds the inverses of the Kii store them in a cell array

if( nargin < 3)
    matOffSet = zeros(ds.chnkrs.npt);
end

opts2 = [];
opts2.adaptive_correction = true;

listKs_inv = cell(1, ds.nGammas);
for i=1:ds.nGammas
    Kii = chunkermat(ds.listGammas(i), kernel, opts2);
    Kii = Kii + matOffSet( (ds.nB(i+1)+1):ds.nB(i+2), (ds.nB(i+1)+1):ds.nB(i+2) );
    listKs_inv{i} = inv(Kii);
end


end