clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point in the following
[~,id]=max(mean(S,1));

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%signal to noise ratio
SNR=1;

%generate noisy data according to given SNR
X=Xs+1/sqrt(SNR)*Noise;

%% MNE: Lambda variation
lambda    = logspace(-1,3,5);
fig_count = 1;   
for k=1:length(lambda)
    fprintf("Figure=%d, Lambda=%d\n", fig_count, lambda(k))
    Shat=MNE(X(:,id),G,lambda(k));
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;
	fig_count = fig_count + 1;
end

%% MNE:Lambda and SNR variation
lambda    = logspace(0,3,4);
SNR       = logspace(-1,1,3);
fig_count = 1;     
for j=1:length(SNR)
    X=Xs+1/sqrt(SNR(j))*Noise;
    for i=1:length(lambda)
        fprintf("Figure=%d, SNR=%d, Lambda=%d\n", fig_count, SNR(j), lambda(i))
        Shat=MNE(X(:,id),G,lambda(i));
        figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;
        fig_count = fig_count + 1;
    end
end
%% MNE: L-curve criterion
SNR    = 1;
X      = Xs+1/sqrt(SNR)*Noise;
lambda = logspace(-10,10,200);
for i=1:length(lambda)
    Shat        = MNE(X(:,id),G,lambda(i));
    err_reco(i) = norm(X(:,id) - G*Shat, 2);
    norm_s(i)   = norm(Shat, 2);
end

figure; loglog(norm_s, err_reco);
title("MNE: L-curve");
ylabel("||X-As||_2");
xlabel("||s||_2");
grid();
set(gca,'fontsize', 24);

[~, idx_err_reco] = min(abs(err_reco - 35.3));
[~, idx_norm_s]   = min(abs(norm_s - 479.6));
fprintf("lambda(||X-As||_2) = %d\n", lambda(idx_err_reco));
fprintf("lambda(||s||_2) = %d\n", lambda(idx_norm_s));
%% MNE: Discrepancy principle
SNR    = 1;
X      = Xs+1/sqrt(SNR)*Noise;
lambda = logspace(-10,10,200);
for i=1:length(lambda)
    Shat        = MNE(X(:,id),G,lambda(i));
    err_reco_2(i) = norm(X(:,id) - G*Shat, 2).^2;
end
norm_n_2 = ones(size(err_reco_2))*norm(1/sqrt(SNR)*Noise, 2).^2;

figure; 
loglog(lambda, err_reco_2, 'DisplayName','||X-As||_2^2'); hold on;
loglog(lambda, norm_n_2, 'DisplayName','||Noise||_2^2'); 
title("MNE: Discrepancy");
ylabel("Power");
xlabel("lambda");
grid();
legend;
set(gca,'fontsize', 24);

fprintf("lambda= %d\n", 1168.5);
%% MNE: Generalized Cross-Validation
%  code to be added in the L-curve iteration
N      = size(G, 1);
SNR    = 1;
X      = Xs+1/sqrt(SNR)*Noise;
lambda = logspace(-10,10,200);
for i=1:length(lambda)
    Shat       = MNE(X(:,id),G,lambda(i));
    err_reco_2 = norm(X(:,id) - G*Shat, 2).^2;
    trace_2    = trace(eye(N) - G*G'*inv(G*G' + lambda(i)*eye(N))).^2;
    GCV(i)     = err_reco_2 / trace_2;
end

figure; loglog(lambda, GCV);
title("MNE: Generalized cross-validation");
ylabel("GCV");
xlabel("lambda");
grid();
set(gca,'fontsize', 24);

[~, idx_GCV] = min(GCV);
fprintf("lambda(GCV) = %d\n", lambda(idx_GCV));

%% MNE: Final choice
SNR    = 1;
X      = Xs+1/sqrt(SNR)*Noise;
lambda = 114.895;
Shat   = MNE(X(:,id),G,lambda);
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;

%% SISSY: Lambda variation
SNR    = 10;
X      = Xs+1/sqrt(SNR)*Noise;
T      = variation_operator(mesh,'face');
alpha = 0.1;
Niter = 60;
lambda = logspace(-2,3,6);
fig_count = 1;
for k=1:length(lambda)
    fprintf("Figure=%d, Lambda=%d\n", fig_count, lambda(k))
    Shat=SISSY(X(:,id),G,T,lambda(k), alpha, Niter);  %  function to be implemented
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;
	fig_count = fig_count + 1;
end

%% SISSY: Alpha variation
SNR    = 10;
X      = Xs+1/sqrt(SNR)*Noise;
T      = variation_operator(mesh,'face');
lambda = 10;
Niter = 60;
alpha = linspace(0, 1, 6);
fig_count = 1;
for k=1:length(alpha)
    fprintf("Figure=%d, Alpha=%d\n", fig_count, alpha(k))
    Shat=SISSY(X(:,id),G,T,lambda, alpha(k), Niter);
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;
	fig_count = fig_count + 1;
end

%% SISSY: L0 restriction criterion
SNR    = 10;
X      = Xs+1/sqrt(SNR)*Noise;
T      = variation_operator(mesh,'face');
alpha = 0.2;
Niter = 60;
lambda = logspace(-5,5,100);
for k=1:length(lambda)
    Shat=SISSY(X(:,id),G,T,lambda(k), alpha, Niter);
    norm_L0_Ts  = sum(T*Shat > 0.01*max(T*Shat));
    norm_L0_s   = sum(Shat > 0.01*max(Shat));
    L0_restriction(k) = norm_L0_Ts + alpha*norm_L0_s;
end

figure; loglog(lambda, L0_restriction);
title("SISSY: L0 restriction criterion");
ylabel("||Ts||_0 + \alpha||s||_0");
xlabel("lambda");
grid();
set(gca,'fontsize', 24);

%% SISSY: Final choice
SNR    = 10;
X      = Xs+1/sqrt(SNR)*Noise;
T      = variation_operator(mesh,'face');
lambda = 5.67243;
alpha  = 0.2;
Niter  = 60;
Shat=SISSY(X(:,id),G,T,lambda, alpha, Niter);
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;