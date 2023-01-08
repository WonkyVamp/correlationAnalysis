function [Priors, Mu, Sigma, Pix] = EM_boundingCov(Data, Priors0, Mu0, Sigma0)
loglik_threshold = 1e-10;
[nbVar, nbData] = size(Data);
nbStates = size(Sigma0,3);
loglik_old = -realmax;
nbStep = 0;

Mu = Mu0;
Sigma = Sigma0;
Priors = Priors0;
while 1
  for i=1:nbStates
    Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
  end
  Pix_tmp = repmat(Priors,[nbData 1]).*Pxi;
  Pix = Pix_tmp ./ repmat(sum(Pix_tmp,2),[1 nbStates]);
  E = sum(Pix);
  for i=1:nbStates
    Priors(i) = E(i) / nbData;
    Mu(:,i) = Data*Pix(:,i) / E(i);
    %Update the covariance matrices
    Data_tmp1 = Data - repmat(Mu(:,i),1,nbData);
    Data_tmp2a = repmat(reshape(Data_tmp1,[nbVar 1 nbData]), [1 nbVar 1]);
    Data_tmp2b = repmat(reshape(Data_tmp1,[1 nbVar nbData]), [nbVar 1 1]);
    Data_tmp2c = repmat(reshape(Pix(:,i),[1 1 nbData]), [nbVar nbVar 1]);
    Sigma(:,:,i) = sum(Data_tmp2a.*Data_tmp2b.*Data_tmp2c, 3) / E(i);
    Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
  end
  for i=1:nbStates
    Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
  end
  F = Pxi*Priors';
  F(find(F<realmin)) = realmin;
  loglik = mean(log(F));
  %Stop the process depending on the increase of the log likelihood 
  if abs((loglik/loglik_old)-1) < loglik_threshold
    break;
  end
  loglik_old = loglik;
  nbStep = nbStep+1;
end

%   for i=1:nbStates
%     Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
%   end
%   for j=1:nbData
%     Pix(j,:) = (Priors.*Pxi(j,:))./(sum(Priors.*Pxi(j,:))+realmin);
%   end
%   E = sum(Pix);
%   for i=1:nbStates
%     Priors(i) = E(i) / nbData;
%     Mu(:,i) = Data*Pix(:,i) / E(i);
%     covtmp = zeros(nbVar,nbVar);
%     for j=1:nbData
%       covtmp = covtmp + (Data(:,j)-Mu(:,i))*(Data(:,j)-Mu(:,i))'.*Pix(j,i);
%     end
%     Sigma(:,:,i) = covtmp / E(i);
%     %% Add a tiny variance to avoid numerical instability
%     Sigma(:,:,i) = Sigma(:,:,i) + 1E-4.*diag(ones(nbVar,1));
%   end
%   for i=1:nbStates
%     %Compute the new probability p(x|i)
%     Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
%   end
%   F = Pxi*Priors';
%   F(find(F<realmin)) = realmin;
%   loglik = mean(log(F));
%   if abs((loglik/loglik_old)-1) < loglik_threshold
%     break;
%   end
%   loglik_old = loglik;
%   nbStep = nbStep+1;
% end


