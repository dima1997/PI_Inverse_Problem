clear;
close all;
clc;

load TP_data;

% ------------ TP1 -----------------------------------
%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point in the following
[~,id]=max(mean(S,1));

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

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

% ------------ end of TP1 ----------------------------


% -----------  TP2 / 3 -------------------------------

%% manual test of lambda
%lambda = logspace(0,3,6);
%for k=1:length(lambda)
%    Shat(:,k)=MNE(X(:,id),G,lambda(k));  % function to be implemented
%    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,k));
%    title(strcat('\lambda=',num2str(lambda(k))));
%end


%% L-curve criterion
% find the L curve for the set of "lambda=logspace(-10,10,20);"

% figure; loglog(n,s);
% n



%% discrepancy principle

%% Generalized Cross-Validation
%  code to be added in the L-curve iteration


%% SISSY
% T=variation_operator(mesh,'face');
%lambda=logspace(-2,3,6); % lambda from 0.01 to 1000
%for k=1:length(lambda)
%[SI,lam]=SISSY(X(:,id),G,T,lambda(k),0.1, MaxIter);  %  function to be implemented
%end
