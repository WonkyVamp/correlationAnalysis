clear all; close all;

%% Load the data for the correct sequences

Data_NN = csvread('../Data/Data_Correct.csv');

addpath('../Utility Functions')

L = size(Data_NN,2);

nDim = 117;

n_seq_corr = size(Data_NN,1)/nDim;

Data = repmat([1:L],1,n_seq_corr);

Data_position=[];
for i=1:n_seq_corr
    Data_position=[Data_position,Data_NN((i-1)*nDim+1:i*nDim,:)];
end
Data=[Data;Data_position];


nbVar = size(Data,1);
nbData = size(Data,2);

nbPC = 4;

[E,v] = eig(cov(Data_position'));
E = fliplr(E);
A = E(:,1:nbPC);
nbVar2 = nbPC+1;
Data_reduced(1,:) = Data(1,:);
Data_reduced(2:nbVar2,:) = A' * Data_position;


nbStates = 5;

[Priors, Mu, Sigma] = EM_init_regularTiming(Data_reduced, nbStates);
[Priors, Mu, Sigma] = EM_boundingCov(Data_reduced, Priors, Mu, Sigma);

disp(['GMM modeling has been completed!',char(10)])


figure('position',[20,120,700,700],'name','GMM Maximum Variance');
for n=1:nbVar2-1
  subplot(nbVar2-1,1,n); hold on;
  plotGMM1(Mu([1,n+1],:), Sigma([1,n+1],[1,n+1],:), [1 0.4 0], 1);
  for j=1:n_seq_corr
     plot(Data_reduced(n+1,(j-1)*L+1:j*L)','color', [0, 0, 0, 0.25], 'LineWidth', 0.25), hold on
  end
  axis([min(Data_reduced(1,:)) max(Data_reduced(1,:)) min(Data_reduced(n+1,:))-0.01 max(Data_reduced(n+1,:))+0.01]);
  set(findobj('type','axes'),'fontsize',10,'box','off')
  xticks(0:20:229)
  yticks(-200:50:200)
  xlabel('Time Frame', 'fontsize',12)
  ylabel('Angle (Degrees)', 'fontsize',12)
end
saveas(gcf, '../Results/PCA_GMM.png')
