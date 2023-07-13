function listKs_inv = computeKii_inv(ds, kernel)
% Builds the inverses of the Kii store them in a cell array


listKs_inv = cell(1, ds.nGammas);
for i=1:ds.nGammas
    Kii = chunkermat(ds.listGammas(i), kernel);
    listKs_inv{i} = inv(Kii);
end


end