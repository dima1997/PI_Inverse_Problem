function [s]=SISSY(x,A,T,lambda,alpha)
    n_iterations=60; rho=1;
    
    z=zeros(size(T,1),size(x,2));
    u=zeros(size(T,1),size(x,2));
    y=zeros(size(T,2),size(x,2));
    v=zeros(size(T,2),size(x,2));
    
    P=sparse(rho*(T.'*T+speye(size(T,2))));
    APi=A/P;
    L=chol(eye(size(A,1))+APi*A.','lower');
    s1=A.'*x;
    
    for k=1:n_iterations
        b=s1+rho*(T.'*(z+u/rho)+y+v/rho);
        s=P\b-APi'*(L'\(L\(APi*b)));
        z=prox_op(T*s-u/rho,'L1',lambda/rho);
        y=prox_op(s-v/rho,'L1',lambda*alpha/rho);
        u=u+rho*(z-T*s);
        v=v+rho*(y-s);
    end
end