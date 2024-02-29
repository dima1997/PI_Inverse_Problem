clear;
close all;
clc;

load TP_data.mat;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('Original source configuration: two source regions','FontSize',18); axis off;

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%signal to noise ratio
SNR=10;

%generate noisy data according to given SNR
X=Xs+1/sqrt(SNR)*Noise;

%visualize data (for a reduced number of sensors whose indices are 
%specified by idx_electrodes)
plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
title('noisy EEG data','FontSize',18);

%% Gibbs sampler

[Shat,Lambdahat] = Gibbs_sampler(Xb(:,id), G, 1/SNR);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:));
title('Reconstructed source configuration','FontSize',18); axis off;
