clear all; close all;


Data_NN_ini = csvread('../Data/Autoencoder_output_correct.csv');

Data_NN_inc_ini = csvread('../Data/Autoencoder_output_incorrect.csv');


scaling_value = 141;
Data_NN = Data_NN_ini*scaling_value;
Data_NN_inc = Data_NN_inc_ini*scaling_value;

addpath('../Utility Functions')

nDim = 4;

L = size(Data_NN,2)/nDim;

n_seq_corr = size(Data_NN,1);

Data=repmat([1:L],1,n_seq_corr);

Data_position=[];
for i=1:n_seq_corr
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN(i,j:nDim:nDim*L)];
    end
    Data_position=[Data_position,temp];
end
Data=[Data;Data_position];

n_seq_inc = size(Data_NN_inc,1);

Data_inc = repmat([1:L],1,n_seq_inc);

Data_inc_position=[];
for i=1:n_seq_corr
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN_inc(i,j:nDim:nDim*L)];
    end
    Data_inc_position=[Data_inc_position,temp];
end
Data_inc = [Data_inc;Data_inc_position];

nbStates = 5;

nbVar = size(Data,1);

[Priors, Mu, Sigma] = EM_init_regularTiming(Data, nbStates);
[Priors, Mu, Sigma] = EM_boundingCov(Data, Priors, Mu, Sigma);


figure('position',[20,120,700,700],'name','GMM Autoencoder');
for n=1:nbVar-1
  subplot(nbVar-1,1,n); hold on;
  plotGMM1(Mu([1,n+1],:), Sigma([1,n+1],[1,n+1],:), [1 0.4 0], 1);
  for j=1:n_seq_corr
     plot(Data(n+1,(j-1)*L+1:j*L)','color', [0, 0, 0, 0.25], 'LineWidth', 0.25), hold on
  end
  axis([min(Data(1,:)) max(Data(1,:)) min(Data(n+1,:))-0.01 max(Data(n+1,:))+0.01]);
  set(findobj('type','axes'),'fontsize',10,'box','off')
  xticks(0:20:229)
  yticks(-200:50:200)
  xlabel('Time Frame', 'fontsize',12)
  ylabel('Angle (Degrees)', 'fontsize',12)
end


for j=1:n_seq_corr
    loglikelihood_corr(j) = loglik(Data(:,(j-1)*L+1:j*L), nbStates, Priors, Mu, Sigma);
end
 
for j=1:n_seq_inc
    loglikelihood_inc(j) = loglik(Data_inc(:,(j-1)*L+1:j*L), nbStates, Priors, Mu, Sigma);
end
 

LL_mean_corr = mean(loglikelihood_corr);

figure, plot(1- abs(loglikelihood_corr - LL_mean_corr)/abs(LL_mean_corr),'ko', 'MarkerFaceColor','g'), hold on,
plot(1- abs(loglikelihood_inc - LL_mean_corr)/abs(LL_mean_corr),'ks', 'MarkerFaceColor','r')
set(findobj('type','axes'),'fontsize',12,'box','off')
axis([0 64 0.4 1.05]);
xlabel('Sequence Number', 'fontsize',14), ylabel('Performance Level', 'fontsize',14);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',12,'location','SW')
saveas(gcf, '../Results/GMM_Performance_Indicators.png')
