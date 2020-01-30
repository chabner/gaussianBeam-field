function [abs_mu] = movmfAbsMu(mu)
% get absolute value from complex mu. gives direction only

abs_mu = real(mu);
abs_mu = abs_mu ./ sqrt(sum(abs(abs_mu).^2,3));

abs_mu_t12 = abs_mu(:,:,1:2);
abs_mu_t3 = abs_mu(:,:,3);

abs_mu_t12(~isfinite(abs_mu_t12)) = 0;
abs_mu_t3(~isfinite(abs_mu_t3)) = 1;

abs_mu = cat(3, abs_mu_t12, abs_mu_t3 );

end

