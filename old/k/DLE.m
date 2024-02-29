function DLE=DLE(idx,idx_est,r_grid)
        r=repmat(r_grid(idx,:),[1,1,size(r_grid(idx_est),1)]);
        r_est=repmat(r_grid(idx_est,:),[1,1,size(r_grid(idx),1)]);
        r_est=permute(r_est,[3 2 1]);
        term_1=mean(min(sqrt(sum((r-r_est).^2,2)),[],3));
        term_2=mean(min(sqrt(sum((r-r_est).^2,2)),[],1));

        DLE=(term_1+term_2)/2;
end