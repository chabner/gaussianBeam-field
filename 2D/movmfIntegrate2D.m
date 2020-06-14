function [I] = movmfIntegrate2D(movmf)
% calculate the integral of movmf
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
%
% OUTPUT
% I: [1,dim(2:end)] integral result

% calculate the real and complex parts


sqrt_mu = sqrt(movmf.mu1.^2 + movmf.mu2.^2);

I_k = 2 * pi * besseli(0,sqrt_mu,1) .* exp(movmf.c + abs(real(sqrt_mu)));

I = sum(movmf.alpha .* I_k,1); % I in size of [1,dim(2:end)]

end

