% SEEMINGLY UNRELATED REGRESSION ESTIMATION
function [Pi, Sigma, Std] = fcn3(X, Y, nx, ny)
XX = (X' * X)^(-1);
XY = X' * Y;
Pi = XX * XY;
MY = Y - X * XX * XY;
Sigma = cov(MY, 1);
Std = Pi ./ ...
  sqrt(repmat(diag(Sigma), 1, nx)' .* repmat(diag(XX), 1, ny));
