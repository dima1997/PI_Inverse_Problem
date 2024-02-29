function [SOut,LambdaOut]=Gibbs_sampler(X,A,sigma2)

%number of iterations
Niter=100;

%get number of sensors and number of dipoles
[N,D]=size(A);

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
        if lambda_i < 0 || lambda_i > 1
            disp("error")
        end
        Q(i,k) = binornd(1,lambda_i, 1, 1);
        if isnan(Q(i,k))
            disp("isnan")
        end
        S(i,k) = Q(i,k) .* normrnd(mu_i, sqrt(sigma_i2)); %((randn(1,1) * sqrt(sigma_i2) + mu_i));
        
        if isinf(S(i,k))
            disp("isinf")
        end
        e = e_i - A(:,i)*S(i,k);
    end
    L = sum(Q(:,k)); 
    lambda = betarnd(1+L,1+D-L);
end 
disp(L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the estimation from Q and S using the MAP criterium
%Q(:,Niter/2:Niter)
cout = zeros(Niter/2,1);
for j = 1:Niter/2
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
disp(size(SOut))
SOut(idx) = (R\(A(:,idx)'*X))/sigma_n2;
disp(size(SOut))
LambdaOut = length(idx)/D; 