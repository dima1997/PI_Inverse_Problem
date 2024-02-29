function [s]=MNE(x,A,lambda)
    s=A'*inv(A*A'+lambda*eye(size(A,1)))*x;
end