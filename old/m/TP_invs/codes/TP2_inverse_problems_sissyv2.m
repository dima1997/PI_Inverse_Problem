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

lambda = 1;
alpha = 0.1;
SNR = 10;
niter = 60;

T = variation_operator(mesh, 'face');

% Generate noisy data according to given SNR
Noise = randn(size(Xs));
X = Xs+1/sqrt(SNR)*Noise;

Shat = SISSY_algo(X(:,id), G, T, lambda, alpha,niter);

figure;
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat, EdgeAlpha=0.2);
view(180,10);
title("SNR:"+SNR+" \lambda:"+lambda+" \alpha:"+alpha,'FontSize',12); 
axis off;

%% lambda
lambda_arr = logspace(-3,4,50);
alpha = 0.2;
SNR = 10;
niter = 60;

%T = variation_operator(mesh, 'face');
Noise = randn(size(Xs));
X = Xs+1/sqrt(SNR)*Noise;

Shat = zeros(D,length(lambda_arr));
for l = 1:length(lambda_arr)
    lambda = lambda_arr(l);
    Shat(:,l) = SISSY_algo(X(:,id), G, T, lambda, alpha,niter);
end

figure; 
k = 1;
for l = 1:length(lambda_arr)
    lambda = lambda_arr(l);
    subplot( 2, 4 ,k);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,l), EdgeAlpha=0.2);
    view(180,10);
    title("\lambda: "+lambda,'FontSize',12); 
    axis off;
    k = k + 1;
end

% Todo continue with the regularization funciton
reg_func = zeros(length(lambda_arr),1);
for l = 1:length(lambda_arr)
    reg_func(l) = zero_norm(T*Shat(:,l)) + alpha*zero_norm(Shat(:,l)); 
end

%% Plot reg func

figure
semilogx(lambda_arr, reg_func)
xlabel("\lambda")

figure
loglog(lambda_arr, reg_func)
grid()
xlabel("\lambda")

%% alpha
lambda = 10;
alpha_arr = linspace(0,1.0,8);
SNR = 10;
niter = 60;

%T = variation_operator(mesh, 'face');
Noise = randn(size(Xs));
X = Xs+1/sqrt(SNR)*Noise;

Shat = zeros(D,length(alpha_arr));
for l = 1:length(alpha_arr)
    alpha = alpha_arr(l);
    Shat(:,l) = SISSY_algo(X(:,id), G, T, lambda, alpha,niter);
end

figure; 
k = 1;
for l = 1:length(alpha_arr)
    alpha = alpha_arr(l);
    subplot( 2, 4 ,k);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,l), EdgeAlpha=0.2);
    view(180,10);
    txt = sprintf("alpha:%3.2f", alpha);
    title(txt,'FontSize',12); 
    axis off;
    k = k + 1;
end
