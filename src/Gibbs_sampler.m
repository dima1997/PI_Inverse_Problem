function [SOut,LambdaOut]=Gibbs_sampler(X, A, sigma_n2, sigma_s2, alpha, beta, Niter)

%get number of sensors and number of dipoles
[~,D]=size(A);

% constants 
nA = sum(A.^2,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% implement the Gibbs sampler here according to the pseudocode
%
% store vectors q and s of each iteration in matrices Q (D x niter) and S
% (D x niter)
% use variables sigma_n2 and sigma_s2 for the variances of noise and
% signals
Q = zeros(D,Niter);
S = zeros(D,Niter);
lambda = betarnd(alpha, beta);
for n = 1:Niter
    fprintf('Gibbs Iter %d\n',n);
    e = X - A * S(:,n);
    for i = 1:D
        e_i = e + A(:,i)*S(i,n);
        sigma_i2 = (sigma_n2 * sigma_s2) / (sigma_n2 + sigma_s2 * nA(i)); 
        mu_i = (sigma_i2 / sigma_n2) * A(:,i)' * e_i;
        nu_i = lambda * (sqrt(sigma_i2) / sqrt(sigma_s2)) * exp((mu_i^2) / (2*sigma_i2));
        if nu_i > 1e10
            lambda_i = 1;
        else
            lambda_i = nu_i / (nu_i + 1 - lambda);
        end
        Q(i,n) = binornd(1,lambda_i, 1, 1);
        S(i,n) = Q(i,n) .* normrnd(mu_i, sqrt(sigma_i2));
        
        e = e_i - A(:,i)*S(i,n);
    end
    L = sum(Q(:,n)); 
    lambda = betarnd(1 + L, 1 + D - L);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the estimation from Q and S using the MAP criterium
cout = zeros(Niter/2,1);
for j = 1:Niter/2
    fprintf('Est. Iter %d \n',j);
    q = Q(:,Niter/2+j); 
    idx = find(q==1); 
    R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
    S_opt = (R\(A(:,idx)'*X))/sigma_n2; 
    cout(j) = - S_opt'*R*S_opt/sigma_n2; 
end

[~, OptIdx] = min(cout); 
q = Q(:,Niter/2+OptIdx);
idx = find(q==1); 
R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
SOut = zeros(D,1); 
SOut(idx) = (R\(A(:,idx)'*X))/sigma_n2;
LambdaOut = length(idx)/D; 