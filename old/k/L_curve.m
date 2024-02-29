function [error,l2_norm]=L_curve(X,G,lambda_range)

    error=zeros(length(lambda_range),1);
    l2_norm=zeros(length(lambda_range),1);

    for k=1:length(lambda_range)
        s_MNE=MNE(X,G,lambda_range(k));
        error(k)=norm(X-G*s_MNE,'fro');
        l2_norm(k)=norm(s_MNE,'fro');
        disp(lambda_range(k));
        disp(l2_norm(k));
    end
end