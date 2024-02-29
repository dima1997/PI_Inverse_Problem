clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

[N, D] = size(G);
%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
%figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
%title('original source configuration: two source regions','FontSize',18); axis off;

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
%plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
%title('noisy EEG data','FontSize',18);

%% recons

lambda_arr = logspace(-10,10,200);

for i = 1:length(lambda_arr)
    lambda = lambda_arr(i);
    shat = MNE_algo(X(:,id),G, lambda);
    error_(i) = norm(X(:,id) - G*shat);
    s_norm(i) = norm(shat);
end
% Shat = Gibbs_sampler(X(:,id),G, 1/SNR);

%% L curve

figure
loglog(s_norm, error_)
grid()
ylabel("||X-As||^2")
xlabel("||s||^2")

[val, idx] = min(abs(s_norm-423));
lambda_arr(idx)
[val, idx] = min(abs(error_-58));
lambda_arr(idx)

%% Discrepancy

figure
loglog(lambda_arr, error_ .^ 2.2)
hold on 
loglog(lambda_arr, ones(size(error_))*norm(Noise)^2, 'r')
grid()
ylabel("||X-As||^2")
xlabel("\lambda")

%% GCV
GCV = error_ .^ 2;
for l = 1:length(lambda_arr)
    lambda = lambda_arr(l);
    GCV(l) = GCV(l) / trace( eye(N) - G * G' / ( G * G' + lambda* eye(N)) ) ^ 2;
end
%% plot
figure
loglog(lambda_arr, GCV)
grid()
ylabel("GCV")
xlabel("\lambda")

%% Plot results 

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:), EdgeAlpha=0.2);
title('Reconstructed source configuration','FontSize',18); axis off;

%% Visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id), EdgeAlpha=0.2);
title('Original source configuration','FontSize',18); axis off;
