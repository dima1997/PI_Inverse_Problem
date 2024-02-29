function [GCV]=generalized_cross_validation(X,G,lambda_range)
    GCV=zeros(length(lambda_range),1);
    for k=1:length(lambda_range)
        s_MNE=MNE(X,G,lambda_range(k));
        GCV(k)=norm(X-G*s_MNE,'fro')^2/(trace(eye(size(G,1))-G*G'*inv(G*G'+lambda_range(k)*eye(size(G,1)))))^2;
    end
end