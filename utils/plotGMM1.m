function plotGMM1(Mu, Sigma, color, display_mode);

nbData = size(Mu,2);
lightcolor = color + [0.6,0.6,0.6];
lightcolor(find(lightcolor>1.0)) = 1.0;

if display_mode==1
  nbDrawingSeg = 40;
  t = linspace(-pi, pi, nbDrawingSeg)';
  for j=1:nbData
    stdev = sqrtm(2.0.*Sigma(:,:,j));
    X = [cos(t) sin(t)] * real(stdev) + repmat(Mu(:,j)',nbDrawingSeg,1);
    patch(X(:,1), X(:,2), lightcolor, 'lineWidth', 2, 'EdgeColor', color);
    plot(Mu(1,:), Mu(2,:), 'x', 'lineWidth', 2, 'color', color);
  end
elseif display_mode==2
  nbDrawingSeg = 40;
  t = linspace(-pi, pi, nbDrawingSeg)';
  for j=1:nbData
    stdev = sqrtm(3.0.*Sigma(:,:,j));
    X = [cos(t) sin(t)] * real(stdev) + repmat(Mu(:,j)',nbDrawingSeg,1);
    patch(X(:,1), X(:,2), lightcolor, 'LineStyle', 'none');
  end
  plot(Mu(1,:), Mu(2,:), '-', 'lineWidth', 3, 'color', color);
elseif display_mode==3
  for j=1:nbData
    ymax(j) = Mu(2,j) + sqrtm(Sigma(1,1,j));  %%% modified
    ymin(j) = Mu(2,j) - sqrtm(Sigma(1,1,j));  %%% modified 
  end
  patch([Mu(1,1:end) Mu(1,end:-1:1)], [ymax(1:end) ymin(end:-1:1)], lightcolor, 'LineStyle', 'none');
  plot(Mu(1,:), Mu(2,:), '-', 'lineWidth', 3, 'color', color); 
end





