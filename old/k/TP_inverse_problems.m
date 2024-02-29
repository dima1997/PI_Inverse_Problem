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
SNR=1;

%generate noisy data according to given SNR
X=Xs+1/sqrt(SNR)*Noise;

%visualize data (for a reduced number of sensors whose indices are 
%specified by idx_electrodes)
plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
title('noisy EEG data','FontSize',18);

%% Gibbs sampler

SNR=logspace(-1,1,3);

for i=1:length(SNR)
    X=Xs+1/sqrt(SNR(i))*Noise;
    subplot(2,2,i)
    [SOut,LambdaOut]=Gibbs_sampler(X(:,id),G); 
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),SOut);
    title('GS with SNR = '+string(SNR(i))); axis off;
end
subplot(2,2,4)
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('Original config'); axis off;

%% MNE algorithm

lambda=1;
s_MNE=MNE(X,G,lambda);
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_MNE(:,id));
title('MNE output: two source regions','FontSize',18); axis off;

%Modifying lambda - SNR=10

lambda=logspace(-1,3,5);

for i=1:length(lambda)
    subplot(3,2,i);
    s_MNE=MNE(X,G,lambda(i));
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_MNE(:,id)); axis off;
    title('MNE output for \lambda = ' + string(lambda(i)))
    hold on
end

subplot(3,2,6)
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id)); axis off;
title('Original config')

%Modifying lambda - Each curve has a fixed SNR

lambda=logspace(0,3,4);
SNR=logspace(-1,1,3);

figure
for j=1:length(SNR)
    for i=1:length(lambda)
        subplot(3,4,i+4*(j-1));
        X=Xs+1/sqrt(SNR(j))*Noise;
        s_MNE=MNE(X,G,lambda(i));
        trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_MNE(:,id)); axis off;
        title('\lambda = ' + string(lambda(i)) + ', SNR = ' + string(SNR(j)))
        hold on
    end
end

% L-curve criterion

SNR=1;
X=Xs+1/sqrt(SNR)*Noise;
lambda=logspace(-10,10,21);

[error,l2_norm]=L_curve(X,G,lambda);

figure
l_curve_range=5:16;
loglog(l2_norm(l_curve_range),error(l_curve_range),'.-')
hold on
text(l2_norm(l_curve_range),error(l_curve_range),string(lambda(l_curve_range)),'FontSize',7)
title('L-curve criterion');
xlabel('Norm of s');
ylabel('Reconstruction error');
dim=[.62 .0 .3 .3];
str='\lambda = '+string(lambda(9));
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% Discrepancy criterion

noise_power=discrepancy_principle(Noise,lambda);

figure
loglog(lambda,noise_power)
hold on
loglog(lambda,error.^2)
title('Discrepancy criterion');
xlabel('\lambda');
ylabel('Reconstruction error');

% Generalized cross validation

GCV=generalized_cross_validation(X,G,lambda);
GCV_range=7:17;
figure
loglog(lambda(GCV_range),GCV(GCV_range))
hold on
title('Generalized cross validation criterion');
xlabel('\lambda');
ylabel('GCV');

%% SISSY algorithm

SNR=10;
X=Xs+1/sqrt(SNR)*Noise;
lambda=1; alpha=0.1;
s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda,alpha);
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_SISSY(:,id));
title('SISSY output: two source regions','FontSize',18); axis off;

%Modifying lambda - SNR=10

lambda=logspace(-2,3,6);

for i=1:length(lambda)
    subplot(3,2,i);
    s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda(i),alpha);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_SISSY(:,id)); axis off;
    title('SISSY output for \lambda = ' + string(lambda(i)))
    hold on
end

%Modifying alpha - lambda=10

alpha=linspace(0,1,5);
lambda=10;

for i=1:length(alpha)
    subplot(3,2,i);
    s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda,alpha(i));
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_SISSY(:,id)); axis off;
    title('\lambda = '+string(lambda)+', \alpha = ' + string(alpha(i)))
    hold on
end
subplot(3,2,6)
trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id)); axis off;
title('Original config')

% Finding lambda using the L0 norm

alpha=0.1;
lambda=logspace(-2,3,6);

for i=1:length(lambda)
        s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda(i),alpha);
        thr=0.01*max(max(s_SISSY));
        s_thr=s_SISSY;
        s_thr(find(abs(s_thr)<abs(thr)))=0;
        Ts_thr=variation_operator(mesh,'face')*s_SISSY;
        thr=0.01*max(max(Ts_thr));
        Ts_thr(find(abs(Ts_thr)<abs(thr)))=0;
        
        constraint(i)=nnz(Ts_thr)+alpha*nnz(s_thr);
end

semilogx(lambda,constraint)
hold on
title('L0-based constraint, \alpha = '+string(alpha)+', SNR = '+string(SNR));
xlabel('\lambda');
ylabel('Constraint value');

%% Qualitative comparison

% Non-correlated noise

SNR=logspace(-1,1,3);
lambda_MNE=1000;
lambda_SISSY=10;
alpha=0.1;

for i=1:length(SNR)
    X=Xs+1/sqrt(SNR(i))*Noise;
    figure
    subplot(1,4,1)
    [SOut,LambdaOut] = Gibbs_sampler(X(:,id),G); 
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),SOut);
    title('GS with SNR = '+string(SNR(i))); axis off;
    hold on
    subplot(1,4,2)
    s_MNE=MNE(X,G,lambda_MNE);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_MNE(:,id)); axis off;
    title('MNE \lambda = '+string(lambda_MNE)+', SNR = '+string(SNR(i)))
    subplot(1,4,3)
    s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda_SISSY,alpha);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_SISSY(:,id)); axis off;
    title('SISSY \lambda = '+string(lambda_SISSY)+', \alpha = '+string(alpha)+', SNR = ' + string(SNR(i)))
    subplot(1,4,4)
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
    title('Original config'); axis off;
end

% Correlated noise

Snoise=zeros(size(S));
Snoise(find(S==0))=randn(size(find(S==0)));
Noise=G*Snoise;
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

for i=1:length(SNR)
    X=Xs+1/sqrt(SNR(i))*Noise;
    figure
    subplot(1,4,1)
    [SOut,LambdaOut] = Gibbs_sampler(X(:,id),G); 
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),SOut);
    title('GS with SNR = '+string(SNR(i))); axis off;
    hold on
    subplot(1,4,2)
    s_MNE=MNE(X,G,lambda_MNE);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_MNE(:,id)); axis off;
    title('MNE \lambda = '+string(lambda_MNE)+', SNR = '+string(SNR(i)))
    subplot(1,4,3)
    s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda_SISSY,alpha);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_SISSY(:,id)); axis off;
    title('SISSY \lambda = '+string(lambda_SISSY)+', \alpha = '+string(alpha)+', SNR = ' + string(SNR(i)))
    subplot(1,4,4)
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
    title('Original config'); axis off;
end

%% Quantitative comparison

% Identifying positions of the original source configuration

r_grid = zeros(size(mesh.f)); 
k=1; 
thr=1;

while(k<=size(mesh.f,1))
    r_grid(k,:)= mesh.v(mesh.f(k,:));
    k=k+1;
end
original=zeros(size(S(:,1)));
inactive=find(S==0);
active=unique(find(abs(S(:,id))>thr));

original(active)=1;
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),original); 
title('Original config'); axis off;

% DLE matrix for each algorithm, SNR=1, 10 sets, spatially correlated noise

DLE_matrix=zeros(10,3);
SNR=1;
lambda_MNE=200;
lambda_SISSY=10;
alpha=0.1;

for it=1:10
    
    Snoise=zeros(size(S));
    Snoise(find(S==0))=randn(size(find(S==0)));
    Noise=G*Snoise;
    Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');
    
    X=Xs+1/sqrt(SNR)*Noise;
    
    [s_GS,LambdaOut]=Gibbs_sampler(X(:,id),G); 
    active_GS=unique(find(abs(s_GS)>thr),size(s_GS,1));

    s_MNE=MNE(X,G,lambda_MNE);
    active_MNE=unique(find(abs(s_MNE(:,id))>thr));
    
    s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda_SISSY,alpha);
    active_SISSY=unique(find(abs(s_SISSY(:,id))>thr));

    DLE_matrix(it,1)=DLE(active,active_GS,r_grid);
    DLE_matrix(it,2)=DLE(active,active_MNE,r_grid);
    DLE_matrix(it,3)=DLE(active,active_SISSY,r_grid);
    
end

boxplot(DLE_matrix,'Labels',{'GS','MNE','SISSY'});
title('DLE boxplot with SNR=1');
grid on;

GS_config=zeros(size(S(:,1)));
GS_config(active_GS)=1;
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),GS_config);
title('GS'); axis off;

MNE_config=zeros(size(S(:,1)));
MNE_config(active_MNE)=1;
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),MNE_config);
title('MNE'); axis off;

SISSY_config=zeros(size(S(:,1)));
SISSY_config(active_SISSY)=1;
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),SISSY_config);
title('SISSY'); axis off;

% DLE for MNE and SISSY, different SNR values, 10 sets

SNR=logspace(-1,1,5);
DLE_matrix=zeros(length(SNR),2);

for it=1:10    
    for j=1:length(SNR)
    
        X=Xs+1/sqrt(SNR(j))*Noise;

        s_MNE=MNE(X,G,lambda_MNE);
        active_MNE=unique(find(abs(s_MNE(:,id))>thr));

        s_SISSY=SISSY(X,G,variation_operator(mesh,'face'),lambda_SISSY,alpha);
        active_SISSY=unique(find(abs(s_SISSY(:,id))>thr));
        active_SISSY(find(active_SISSY==0))=size(active_SISSY,1);
        
        DLE_matrix(j,1)=DLE_matrix(j,1)+DLE(active,active_MNE,r_grid)/10;
        DLE_matrix(j,2)=DLE_matrix(j,2)+DLE(active,active_SISSY,r_grid)/10;
    
    end
       
end

semilogx(SNR,DLE_matrix(:,1));
hold on
semilogx(SNR,DLE_matrix(:,2));
xlabel('Signal to Noise Ratio');
ylabel('Mean DLE');
legend('MNE','SISSY');
title('Mean DLE as a function of the SNR');
grid on;