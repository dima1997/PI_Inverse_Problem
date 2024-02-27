function [SOut,LambdaOut]=Gibbs_sampler(X,A, sigma_n2, sigma_s2, alpha, beta, Niter)

%number of iterations
##Niter=100;

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
Q = zeros(D, Niter);
S = zeros(D, Niter);
##sigma_n2 = 1;
##sigma_s2 = 1;
##alpha    = 1;
##beta     = 1;
lambda   = betarnd(alpha, beta);
for n = 1:Niter
  fprintf('Iteration %d\n', n);
  e = X - A*S(:,n);
  for i = 1:D
    e_i      = e + A(:,i)*S(i,n);
    sigma_i2 = (sigma_n2*sigma_s2)/(sigma_n2 + sigma_s2*nA(i));
    mu_i     = (sigma_i2/sigma_n2)*A(:,i)'*e_i;
    nu_i     = lambda*(sqrt(sigma_i2)/sqrt(sigma_s2))*exp((mu_i^2)/(2*sigma_i2));
    lambda_i = nu_i / (nu_i + 1 + lambda);
    Q(i,n)   = binornd(1, lambda_i);
    S(i,n)   = Q(i,n) * normrnd(mu_i, sqrt(sigma_i2));
    e        = e_i - A(:,i)*S(i,n);
  endfor
  L      = sum(Q(:,n));
  lambda = betarnd(alpha + L, beta + D - L);
endfor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the estimation from Q and S using the MAP criterium
%Q(:,Niter/2:Niter)

for j = 1:Niter/2
    q = Q(:,Niter/2+j);
    idx = find(q==1);
    R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2;
    S_opt = (R\(A(:,idx)'*X))/sigma_n2;
    cout(j) = - S_opt'*R*S_opt/sigma_n2;
end

[Val, OptIdx] = min(cout);
q = Q(:,Niter/2+OptIdx);
idx = find(q==1);
R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2;
SOut = zeros(D,1);
SOut(idx) = (R\(A(:,idx)'*X))/sigma_n2;
LambdaOut = length(idx)/D;



