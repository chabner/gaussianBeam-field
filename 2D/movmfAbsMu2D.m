function [abs_mu1,abs_mu2,kappa] = movmfAbsMu2D(movmf)
% get absolute value from complex mu. gives direction only

abs_mu1 = real(movmf.mu1);
abs_mu2 = real(movmf.mu2);

kappa = sqrt(abs_mu1.^2 + abs_mu2.^2);

abs_mu1 = abs_mu1 ./ kappa;
abs_mu2 = abs_mu2 ./ kappa;

abs_mu1(~isfinite(abs_mu1)) = 0;
abs_mu2(~isfinite(abs_mu2)) = 1;

end


