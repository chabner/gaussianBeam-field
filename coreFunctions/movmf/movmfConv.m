function [conv_vmf]=movmfConv(apertureVmf,scatteringMovmf)
% convolve aperture with scattering function

% INPUT
% apertureVmf: complex vmf aperture.
% scatteringMovmf: vmf mixture of the scattering function, mu only in
% directions of [0,0,1] or [0,0,-1]
%
% OUTPUT
% conv_vmf: convolution result.

[mu_r1,mu_r2,mu_r3] = movmfAbsMu(apertureVmf);

% build rotation matrix that rotates mu_r to the north pole
Rs = (1 ./ (1 + mu_r3));

R_11 = 1 - (mu_r1.^2) .* Rs;
R_12 = -mu_r1 .* mu_r2 .* Rs;
R_13 = -mu_r1;

R_21 = R_12;
R_22 = 1 - (mu_r2.^2) .* Rs;
R_23 = -mu_r2;

R_31 = mu_r1;
R_32 = mu_r2;
R_33 = 1 - (mu_r1.^2 + mu_r2.^2) .* Rs;

% aperture to gaussian
rotatedMu_1 = R_11 .* apertureVmf.mu1 + R_12 .* apertureVmf.mu2 + R_13 .* apertureVmf.mu3;
rotatedMu_2 = R_21 .* apertureVmf.mu1 + R_22 .* apertureVmf.mu2 + R_23 .* apertureVmf.mu3;
rotatedMu_3 = R_31 .* apertureVmf.mu1 + R_32 .* apertureVmf.mu2 + R_33 .* apertureVmf.mu3;

A1 = rotatedMu_3;

B1_1 = rotatedMu_1;
B1_2 = rotatedMu_2;
c1 = rotatedMu_3 + apertureVmf.c;

% convert to gaussian
A2 = scatteringMovmf.mu3;
c2 = scatteringMovmf.mu3 + scatteringMovmf.c;

% convolution + gaussian back to vmf
sumA1A2 = A1 + A2;
invSumA1A2 =  1 ./ sumA1A2;

conv_mu3 = A2 .* (invSumA1A2 .* A1);

conv_mu1 = (B1_1 .* invSumA1A2) .* A2;
conv_mu2 = (B1_2 .* invSumA1A2) .* A2;

conv_c = c1 + c2 + (B1_1.^2 + B1_2.^2) .* invSumA1A2 / 2 + log((2*pi) ./ sumA1A2) - conv_mu3;

% rotate back
conv_vmf.mu1 = R_11 .* conv_mu1 + R_21 .* conv_mu2 + R_31 .* conv_mu3;
conv_vmf.mu2 = R_12 .* conv_mu1 + R_22 .* conv_mu2 + R_32 .* conv_mu3;
conv_vmf.mu3 = R_13 .* conv_mu1 + R_23 .* conv_mu2 + R_33 .* conv_mu3;

% build the vmf
conv_vmf.alpha = repmat(scatteringMovmf.alpha,1,apertureVmf.dim(2));
conv_vmf.c = conv_c;

conv_vmf.dim = size(conv_vmf.alpha);

% normlize to ideal result
% take absolute, maximal direction is mu_r
[w_max_1,w_max_2,w_max_3] = movmfAbsMu(conv_vmf);

% maximal value in maximal direction
log_estimated_conv_max = conv_vmf.mu1 .* w_max_1 + ...
                         conv_vmf.mu2 .* w_max_2 + ...
                         conv_vmf.mu3 .* w_max_3 + conv_vmf.c;
                     
C = sqrt((apertureVmf.mu1 + abs(scatteringMovmf.mu3) .* w_max_1).^2 + ...
         (apertureVmf.mu2 + abs(scatteringMovmf.mu3) .* w_max_2).^2 + ...
         (apertureVmf.mu3 + abs(scatteringMovmf.mu3) .* w_max_3).^2);

c = apertureVmf.c + scatteringMovmf.c;
log_exact_conv_max = c+C + log(2*pi) - log(C);

% fix the maximum where it matters
conv_vmf.c = conv_vmf.c - log_estimated_conv_max + log_exact_conv_max;

end

