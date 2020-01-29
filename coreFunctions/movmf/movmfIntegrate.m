function [I] = movmfIntegrate(movmf)
% calculate the integral of movmf
%
% INPUT
% movmf: mixture of von Mises-Fisher structre
%
% OUTPUT
% I: [N,1] integral result

% calculate the real and complex parts

sqrt_mu = sqrt(sum(movmf.mu.^2,3));

I_k = 2 * pi * exp(movmf.c + sqrt_mu) ./ (sqrt_mu);
% I_k(sqrt_mu == 0) = exp(movmf.c(sqrt_mu == 0)) .* 4 * pi;

I = sum(movmf.alpha .* I_k,2); % I in size of [N,1]

end

