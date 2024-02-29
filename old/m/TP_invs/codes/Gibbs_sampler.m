function [SOut,LambdaOut]=Gibbs_sampler(X,A,sigma2)

%number of iterations
Niter=100;
Nsearch = Niter/2;

%get number of sensors and number of dipoles
[N,D]=size(A);

% constants 
nA = sum(A.^2,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implement the Gibbs sampler here according to the pseudocode
% store vectors q and s of each iteration in matrices Q and S (D x niter)
% use variables sigma_n2 and sigma_s2 for the variances of noise and s
Q = zeros(D,Niter);
S = zeros(D,Niter);
lambda = betarnd(1,1);
sigma_n2 = sigma2;
sigma_s2 = sigma2;
for k = 1:Niter
    fprintf('Iteration: %d \n',k);
    e = X - A * S(:,k);
    for i = 1:D
        e_i = e + A(:,i)*S(i,k);
        sigma_i2 = sigma_n2 * sigma_s2 / (sigma_n2 + sigma_s2 * nA(i) ); %sum(A(:,i)^2) nA(i)
        mu_i = (sigma_i2 / sigma_n2) * A(:,i)' * e_i;
        v_i = lambda * sqrt(sigma_i2 / sigma_s2) * exp(mu_i^2 / (2*sigma_i2));
        if v_i > 1e10
            lambda_i = 1;
        else
            lambda_i = v_i / (v_i + 1 - lambda);
        end
        Q(i,k) = binornd(1,lambda_i, 1, 1);
        S(i,k) = Q(i,k) .* normrnd(mu_i, sqrt(sigma_i2)); %((randn(1,1) * sqrt(sigma_i2) + mu_i));        
        e = e_i - A(:,i)*S(i,k);
    end
    L = sum(Q(:,k)); 
    lambda = betarnd(1+L,1+D-L);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the estimation from Q and S using the MAP criterium
%Q(:,Niter/2:Niter)
cout = zeros(Nsearch,1);
for j = 1:Nsearch
    q = Q(:,Niter-Nsearch+j); 
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