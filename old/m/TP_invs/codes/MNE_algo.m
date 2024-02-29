function s = MNE_algo(X, A, lambda)
    [N, ~] = size(A);
    s = A' * ((A * A' + lambda * eye(N)) \ X);
end 