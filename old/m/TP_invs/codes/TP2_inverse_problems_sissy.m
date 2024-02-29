clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');





%% Noise addition and  Reconstruction
%signal to noise ratio
%SNR_arr = [0.1, 1, 10];
%lambda_arr = [0.1, 1, 10, 100, 1000];

% For the Lcurve
SNR_arr = [0.1, 1, 10];
lambda_arr = logspace(-1,4,100);

[~, D] = size(G);

Shat = zeros(D, length(lambda_arr), length(SNR_arr));

for i = 1:length(SNR_arr)
    SNR = SNR_arr(i);

    %generate noisy data according to given SNR
    X=Xs+1/sqrt(SNR)*Noise;        
    for l = 1:length(lambda_arr)
        lambda = lambda_arr(l);
        
        Shat(:,l,i) = MNE_algo(X(:,id), G, lambda);
    end
end
%% Plot all results 

figure; 
k = 1;
for i = 1:length(SNR_arr)
    SNR = SNR_arr(i);
    for l = 1:length(lambda_arr)
        lambda = lambda_arr(l);
    
        
        subplot(length(SNR_arr), length(lambda_arr),k);
        trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,l,i), EdgeAlpha=0.2);
        view(180,10);
        title("SNR:"+SNR+" lambda:"+lambda,'FontSize',12); 
        axis off;
        k = k + 1;
    end
end


%% Visualize original source distribution

l_idx = 2;
i_idx = 2;
figure;
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,l_idx,i_idx), EdgeAlpha=0.2);
view(180,10);
title("SNR:"+SNR_arr(i_idx)+" lambda:"+lambda_arr(l_idx),'FontSize',12); 
axis off;

figure; 
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id), EdgeAlpha=0.2);
view(180,10);
title('Original source configuration','FontSize',18); axis off;

%% Visualize one reult
lambda = 1;
alpha = 0.1;
SNR = 10;
niter = 60;

T = variation_operator(mesh, 'face');


%generate noisy data according to given SNR
Noise=randn(size(Xs));
X=Xs+1/sqrt(SNR)*Noise;   
Shat = SISSY_algo(X(:,id), G, T, lambda, alpha,niter);

figure;
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat, EdgeAlpha=0.2);
view(180,10);
title("SNR:"+SNR_arr(i_idx)+" lambda:"+lambda_arr(l_idx),'FontSize',12); 
axis off;
        