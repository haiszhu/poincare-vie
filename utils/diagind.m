function i = diagind(A) % return indices of diagonal of square matrix
N = size(A,1); i = sub2ind(size(A), 1:N, 1:N);
end
