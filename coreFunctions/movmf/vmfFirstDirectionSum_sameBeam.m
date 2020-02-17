function [w0,w0p] = vmfFirstDirectionSum_sameBeam(movmf,n)

dim = 3;

% normalize
movmf = movmfAbs(movmf);

% square
movmf.c = 2 * movmf.c(:,n);
movmf.mu1 = 2 * movmf.mu1(:,n);
movmf.mu2 = 2 * movmf.mu2(:,n);
movmf.mu3 = 2 * movmf.mu3(:,n);
movmf.alpha = movmf.alpha(:,n).^2;

[mu_r1, mu_r2, mu_r3, kappa] = movmfAbsMu(movmf);

log_C = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(kappa);
log_C(kappa == 0) = - log(4*pi);

c_1_2 = real(movmf.c);
alpha_1_2 = abs(movmf.alpha);

alpha = alpha_1_2 .* exp(c_1_2 - log_C);

movmf.alpha = alpha / sum(alpha);
movmf.c = log_C;

% sample a direction
smpNum = datasample(1:1:movmf.dim(1),1,'Weights',gather(movmf.alpha));
w0 = vsamp(gather([mu_r1(smpNum);mu_r2(smpNum);mu_r3(smpNum)]), gather(kappa(smpNum)), 1);

% calculate probability
w0p = sqrt(movmfPdf(movmf,w0));

w0 = w0(:);

end

