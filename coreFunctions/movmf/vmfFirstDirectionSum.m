function [w0,w0p] = vmfFirstDirectionSum(movmf)

dim = 3;

% normalize
movmf = movmfAbs(movmf);

% square
movmf.c = 2 * movmf.c;
movmf.mu1 = 2 * movmf.mu1;
movmf.mu2 = 2 * movmf.mu2;
movmf.mu3 = 2 * movmf.mu3;
movmf.alpha = movmf.alpha.^2;

[mu_r1, mu_r2, mu_r3, kappa] = movmfAbsMu(movmf);

log_C = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(kappa);
log_C(kappa == 0) = - log(4*pi);

c_1_2 = real(movmf.c);
alpha_1_2 = abs(movmf.alpha);

alpha = alpha_1_2 .* exp(c_1_2 - log_C);

movmf.alpha = alpha ./ sum(alpha,2);
movmf.c = log_C;

n = randi(movmf.dim(2));

% sample a direction
smpNum = datasample(1:1:movmf.dim(1),1,'Weights',gather(movmf.alpha(:,n)));
w0 = vsamp(gather([mu_r1(smpNum,n);mu_r2(smpNum,n);mu_r3(smpNum,n)]), gather(kappa(smpNum,n)), 1);

% calculate probability
w0p = sqrt(mean(movmfPdf(movmf,w0)));

w0 = w0(:);

end

