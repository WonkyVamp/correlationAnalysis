function [Priors, Mu, Sigma] = EM_init_regularTiming(Data, nbStates)

TimingSep = linspace(min(Data(1,:)), max(Data(1,:)), nbStates+1);

for i=1:nbStates
  idtmp = find( Data(1,:)>=TimingSep(i) & Data(1,:)<TimingSep(i+1));
  Priors(i) = length(idtmp);
  Mu(:,i) = mean(Data(:,idtmp)');
  Sigma(:,:,i) = cov(Data(:,idtmp)');
end
Priors = Priors ./ sum(Priors);


