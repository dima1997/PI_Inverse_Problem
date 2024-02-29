clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('original source configuration: two source regions','FontSize',18); axis off;

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

% Shat = Gibbs_sampler(X(:,id),G, 1/SNR);

%% Plot results 

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:), EdgeAlpha=0.2);
title('Reconstructed source configuration','FontSize',18); axis off;

%% Visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id), EdgeAlpha=0.2);
title('Original source configuration','FontSize',18); axis off;


