%% TP1
clear;
close all;
clc;

load TP_data;

%% Data
%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point in the following
[~,id]=max(mean(S,1));

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%% Original source
%visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('Original source','FontSize',18); axis off;

%% Gibbs sampler
%signal to noise ratios
SNRs = [0.1, 1, 10]; 
for i = 1:length(SNRs)
    %signal to noise ratio
    SNR = SNRs(i);
    fprintf("SNR = %d\n", SNR);
    %generate noisy data according to given SNR
    X=Xs+1/sqrt(SNR)*Noise;
    sigma_n2 = 1/SNR; 
    sigma_s2 = 1/SNR; 
    alpha    = 1; 
    beta     = 1; 
    NIter    = 100;
    [Shat,Lambdahat] = Gibbs_sampler(X(:,id),G,sigma_n2, sigma_s2, alpha, beta, NIter);
    %visualize data (for a reduced number of sensors whose indices are
    %specified by idx_electrodes)
    plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
    title(sprintf('EEG data (SNR=%d)', SNR),'FontSize',18);
    %visualize reconstructed source distribution
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:));
    title(sprintf('Gibbs sampling (SNR=%d)', SNR),'FontSize',18); axis off;
end