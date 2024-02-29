function  s = SISSY_algo(X, A, T,lambda, alpha, niter)
    rho = 1;
    [~,D] = size(A);
    %s = zeros(D,1);
    y = zeros(D,1);
    z = zeros(size(T,1),1);
    u = zeros(size(T,1),1);
    v = zeros(D,1);
    
    P=sparse(rho*(T.'*T+speye(size(T,2))));
    APi=A/P;
    L=chol(eye(size(A,1))+APi*A.',"lower");
    s1=A.'*X;
    
    %P=sparse(rho*(T'*T+speye(size(T,2))));
    %APi = A/P;
    %L=chol(eye(size(A,1))+APi*A',"lower");
    %s = A'*X;
    
    for i = 1:niter
        %s = (A' * A + rho * (T'*T + eye(D))) \ (A'*X + rho * T' * z + T' * u + rho * y + v);

        b=s1+rho*(T.'*(z+u/rho)+y+v/rho);
        s=P\b-APi'*(L'\(L\(APi*b)));
        
        %b = s + rho*(T' * (z + u/rho) + y + v/rho);
        %s = P\b-APi'*(L'\(L\(APi*b)));
        
        
        z = prox_op( T*s - u / rho ,'L1' , (lambda)/rho);
        y = prox_op( s - v / rho,'L1' , (lambda*alpha)/rho);
        u = u + rho*(z - T*s);
        v = v + rho*(y - s);
    end
end