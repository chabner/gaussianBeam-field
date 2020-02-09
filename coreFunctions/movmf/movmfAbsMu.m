function [abs_mu1,abs_mu2,abs_mu3,kappa] = movmfAbsMu(movmf)
% get absolute value from complex mu. gives direction only

abs_mu1 = real(movmf.mu1);
abs_mu2 = real(movmf.mu2);
abs_mu3 = real(movmf.mu3);

kappa = sqrt(abs_mu1.^2 + abs_mu2.^2 + abs_mu3.^2);

abs_mu1 = abs_mu1 ./ kappa;
abs_mu2 = abs_mu2 ./ kappa;
abs_mu3 = abs_mu3 ./ kappa;

abs_mu1(~isfinite(abs_mu1)) = 0;
abs_mu2(~isfinite(abs_mu2)) = 0;
abs_mu3(~isfinite(abs_mu3)) = 1;

end


