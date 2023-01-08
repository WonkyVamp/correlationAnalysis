function prob = gaussPDF(Data, Mu, Sigma)

[nbVar,nbData] = size(Data);

Data = Data' - repmat(Mu',nbData,1);
prob = sum((Data*inv(Sigma)).*Data, 2);
prob = exp(-0.5*prob) / sqrt((2*pi)^nbVar * (abs(det(Sigma))+realmin));
