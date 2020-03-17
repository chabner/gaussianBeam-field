function [I] = movmfIntegrate(movmf)
% calculate the integral of movmf
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
%
% OUTPUT
% I: [1,dim(2:end)] integral result

% calculate the real and complex parts


sqrt_mu = sqrt(movmf.mu1.^2 + movmf.mu2.^2 + movmf.mu3.^2);

I_k = 2 * pi * exp(movmf.c + sqrt_mu) ./ (sqrt_mu);
% I_k(sqrt_mu == 0) = exp(movmf.c(sqrt_mu == 0)) .* 4 * pi;

I = sum(movmf.alpha .* I_k,1); % I in size of [1,dim(2:end)]

end

