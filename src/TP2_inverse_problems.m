clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point in the following
[~,id]=max(mean(S,1));

%visualize original source distribution
% figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
% title('original source configuration: two source regions','FontSize',18); axis off;

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%signal to noise ratio
SNR=1;

%generate noisy data according to given SNR
X=Xs+1/sqrt(SNR)*Noise;

%visualize data (for a reduced number of sensors whose indices are
%specified by idx_electrodes)
% plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
% title('noisy EEG data','FontSize',18);

%% Lambda variation
lambda    = logspace(-1,3,5);
fig_count = 1;   
for k=1:length(lambda)
    fprintf("Figure=%d, Lambda=%d\n", fig_count, lambda(k))
    Shat=MNE(X(:,id),G,lambda(k));
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat); axis off;
	fig_count = fig_count + 1;
end

%% Lambda and SNR variation
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
%% L-curve criterion
SNR    = 1;
X      = Xs+1/sqrt(SNR)*Noise;
lambda = logspace(-10,10,200);
for i=1:length(lambda)
    Shat        = MNE(X(:,id),G,lambda(i));
    err_reco(i) = norm(X(:,id) - G*Shat, 'fro');
    norm_s(i)   = norm(Shat, 'fro');
end

figure; loglog(norm_s, err_reco);
title("MNE: L-curve");
ylabel("||X-As||");
xlabel("||s||");
grid();
set(gca,'fontsize', 24);

[~, idx_err_reco] = min(abs(err_reco - 2.34));
[~, idx_norm_s]   = min(abs(norm_s - 283.3));
fprintf("lambda(||X-As||) = %d\n", lambda(idx_err_reco))
fprintf("lambda(||s||) = %d\n", lambda(idx_norm_s))
%% discrepancy principle
SNR    = 1;
X      = Xs+1/sqrt(SNR)*Noise;
lambda = logspace(-10,10,200);
for i=1:length(lambda)
    Shat        = MNE(X(:,id),G,lambda(i));
    err_reco(i) = norm(X(:,id) - G*Shat, 'fro');
end
norm_n = ones(size(err_reco))*norm(Noise, 'fro');

figure; 
loglog(lambda, err_reco.^2.75, 'DisplayName','||X-As||^2'); hold on;
loglog(lambda, norm_n.^2, 'DisplayName','||Noise||^2'); 
title("MNE: Discrepancy");
ylabel("Power");
xlabel("lambda");
grid();
legend;
set(gca,'fontsize', 24);

fprintf("lambda= %d\n", 11758.5)
%% Generalized Cross-Validation
%  code to be added in the L-curve iteration


%% SISSY
T=variation_operator(mesh,'face');
lambda=logspace(-2,3,6); % lambda from 0.01 to 1000
for k=1:length(lambda)
[SI,lam]=SISSY(X(:,id),G,T,lambda(k),0.1, MaxIter);  %  function to be implemented
end
