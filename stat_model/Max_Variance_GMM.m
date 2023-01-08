
clear all; close all;

%% Load the data for the correct sequences

Data_NN = csvread('../Data/Data_Correct.csv');

addpath('../Utility Functions')

n_dim_reduce = 4;

L = size(Data_NN,2);

nDim = 117;

n_seq_corr = size(Data_NN,1)/nDim;

stdev = std(Data_NN(1:117,:)',1);
stdev_sorted = sort(stdev,'descend');
for i=1:n_dim_reduce
    ind(i) = find(stdev == stdev_sorted(i));
end

Data=repmat([1:L],1,n_seq_corr);

Data_position=[];
for i = 1:n_seq_corr
    repetition_temp = Data_NN((i-1)*nDim+1:i*nDim,:);
    % Extract the 4 dimensions with largest variatioin
    Data_reduced_temp = repetition_temp(ind,:);
    Data_position=[Data_position,Data_reduced_temp];
end
Data=[Data;Data_position];

nbStates = 5;

nbVar = size(Data,1);

[Priors, Mu, Sigma] = EM_init_regularTiming(Data, nbStates);
[Priors, Mu, Sigma] = EM_boundingCov(Data, Priors, Mu, Sigma);

disp(['GMM modeling has been completed!',char(10)])


figure('position',[20,120,700,700],'name','GMM PCA');
for n=1:nbVar-1
  subplot(nbVar-1,1,n); hold on;
  plotGMM_Calinon(Mu([1,n+1],:), Sigma([1,n+1],[1,n+1],:), [1 0.4 0], 1);
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
saveas(gcf, '../Results/Max_Variance_GMM.png')

