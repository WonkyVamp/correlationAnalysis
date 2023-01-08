function [ ll ] = loglik( Data, nbStates, Priors, Mu, Sigma)

for i=1:nbStates
    Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
end

F = Pxi*Priors';
F(find(F<realmin)) = realmin;
ll = sum(log(F))/size(Data,2);

end

