function [SOut,LambdaOut]=Gibbs_sampler(X,A)

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

Q=zeros(D,Niter);
S=zeros(D,Niter);
q=zeros(D,1);
s=zeros(D,1);

%lambda=0.0101; % nnz(S)/(size(S,1)*size(S,2))
sigma_n2=3e3;
sigma_s2=20;
lambda=betarnd(1.01,2e3); 

for it=1:Niter
    e=X-A*s;
    for i=1:D

       ei=e+A(:,i)*s(i); 
       sigma_i=sigma_n2*sigma_s2/(sigma_n2+sigma_s2*norm(A(:,i))^2);
       mu_i=sigma_i/sigma_n2*A(:,i)'*ei;
       v_i=lambda*sqrt(sigma_i/sigma_s2)*exp(mu_i^2/(2*sigma_i)); 
       lambda_i=v_i/(v_i+1-lambda);

       if isnan(lambda_i)
           qi=1;
       else
           qi=binornd(1,lambda_i);
       end

       if qi==0
           si=0;
       else
           si=mu_i+sqrt(sigma_i)*randn();
       end
       
       e=ei-A(:,i)*si; 

       Q(i,it)=qi;
       S(i,it)=si;
    end
    
    q=Q(:,it);
    s=S(:,it);
    
    L=sum(q);
    lambda=betarnd(1.01+L,2e3+D-L); 
    %lambda=L/D;

end
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



