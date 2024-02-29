function [SOut]=MNE(x,A,lambda)

[N,~]=size(A);

SOut = A' * inv(A * A' + lambda * eye(N)) * x;
