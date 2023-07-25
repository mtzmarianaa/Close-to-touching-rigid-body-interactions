function mat = evaluateRInterpolated(xDist, matInterpolantChunk)
% Given a nxnxk matrix with coefficients evaluate the interpolation at
% xDist (number)

n = size(matInterpolantChunk, 1);
mat = zeros( n,n ); % Initialize matrix
k = size(matInterpolantChunk, 3);

for row=1:n
    for col=1:n
        coefs = reshape(matInterpolantChunk(row, col, :), k, 1);
        mat(row, col) = lege.exev(xDist, coefs);
    end
end


end