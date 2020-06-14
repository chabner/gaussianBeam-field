function [w0,w0p] = vmfFirstDirectionSum_sameBeam(movmf,n)

dim = 3;

% normalize
movmf = movmfAbs(movmf);
aperturesNum = prod(movmf.dim(2:end));
totalSmpNum = numel(n);

N = n(:) + aperturesNum * (0:1:(totalSmpNum-1)).';

c = movmf.c(:,N);
mu1 = movmf.mu1(:,N);
mu2 = movmf.mu2(:,N);
mu3 = movmf.mu3(:,N);
alpha = movmf.alpha(:,N);

% square
movmf.c = reshape(permute(c,[1,3,2]) + permute(c,[3,1,2]),[],size(c,2));
movmf.mu1 = reshape(permute(mu1,[1,3,2]) + permute(mu1,[3,1,2]),[],size(mu1,2));
movmf.mu2 = reshape(permute(mu2,[1,3,2]) + permute(mu2,[3,1,2]),[],size(mu2,2));
movmf.mu3 = reshape(permute(mu3,[1,3,2]) + permute(mu3,[3,1,2]),[],size(mu3,2));
movmf.alpha = reshape(permute(alpha,[1,3,2]) .* permute(alpha,[3,1,2]),[],size(alpha,2));

[mu_r1, mu_r2, mu_r3, kappa] = movmfAbsMu(movmf);

log_C = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim,kappa);
log_C(kappa == 0) = - log(4*pi);

c_1_2 = real(movmf.c);
alpha_1_2 = abs(movmf.alpha);

alpha = alpha_1_2 .* exp(c_1_2 - log_C);

movmf.alpha = alpha ./ sum(alpha,1);
movmf.c = log_C;

w0 = zeros(3,totalSmpNum);

% sample a direction
for currSmp = 1:1:totalSmpNum
    smpNum = datasample(1:1:(size(movmf.alpha,1)),1,'Weights',movmf.alpha(:,currSmp));
    w0(:,currSmp) = vsamp([mu_r1(smpNum,currSmp);mu_r2(smpNum,currSmp);mu_r3(smpNum,currSmp)], kappa(smpNum,currSmp), 1);
end

% calculate probability
muTimesX = movmf.mu1.*w0(1,:) + movmf.mu2.*w0(2,:) + movmf.mu3.*w0(3,:);
w0p = sqrt(sum(movmf.alpha .* exp(muTimesX + movmf.c) , 1));

w0 = reshape(w0,[dim,1,1,1,1,totalSmpNum]);
w0p = reshape(w0p,[1,1,1,1,1,totalSmpNum]);
end

