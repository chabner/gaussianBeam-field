function [conv_vmf]=movmfConv2D(apertureVmfL,scatteringMovmf,dirv,dirDim)
% convolve aperture with scattering function

% INPUT
% apertureVmf: complex vmf aperture.
% scatteringMovmf: vmf mixture of the scattering function, mu only in
% directions of [0,1] or [0,-1]
%
% OUTPUT
% conv_vmf: convolution result.
perm = 1:1:4; perm(2) = dirDim; perm(dirDim) = 2;

mu_v1 = permute(dirv(1,:),perm); 
mu_v2 = permute(dirv(2,:),perm); 

[mu_r1,mu_r2] = movmfAbsMu2D(apertureVmfL);

% build rotation matrix that rotates dirv to [0,1]

R_11 = mu_v2;
R_12 = -mu_v1;

R_21 = mu_v1;
R_22 = mu_v2;

% multiple the aperture with the transpose of the rotation matrix 
rotatedMu_1 = R_11 .* apertureVmfL.mu1 + R_21 .* apertureVmfL.mu2;
rotatedMu_2 = R_12 .* apertureVmfL.mu1 + R_22 .* apertureVmfL.mu2;

kappaG = abs(scatteringMovmf.mu2);

conv_const = kappaG ./ (kappaG + rotatedMu_2);
% convSign = sign(real(conv_const.*(rotatedMu_1.*mu_r1_rotated + rotatedMu_2.*mu_r2_rotated + rotatedMu_3.*mu_r3_rotated)));
% convSign(convSign == 0) = 1;

mu1 = conv_const .* rotatedMu_1;
mu2 = conv_const .* rotatedMu_2;

conv_vmf_sol1.mu1 = R_11 .* mu1 + R_12 .* mu2;
conv_vmf_sol1.mu2 = R_21 .* mu1 + R_22 .* mu2;

conv_vmf_sol1.dim = size(conv_vmf_sol1.mu2);
conv_vmf_sol1.alpha = repmat(scatteringMovmf.alpha,[1,conv_vmf_sol1.dim(2:end)]);
conv_vmf_sol1.dim(end+1:4) = 1;

% normlize to ideal result
% take absolute, maximal direction is mu_r
[w_max_1,w_max_2] = movmfAbsMu2D(conv_vmf_sol1);

% maximal value in maximal direction
log_estimated_conv_max = conv_vmf_sol1.mu1 .* w_max_1 + ...
                         conv_vmf_sol1.mu2 .* w_max_2;
                     
C = sqrt((apertureVmfL.mu1 + kappaG .* w_max_1).^2 + ...
         (apertureVmfL.mu2 + kappaG .* w_max_2).^2);

c = apertureVmfL.c + scatteringMovmf.c;
log_exact_conv_max = c+log(besseli(0,C,1)) + abs(real(C)) + log(2*pi);

% fix the maximum where it matters
conv_vmf_sol1.c = - log_estimated_conv_max + log_exact_conv_max;


sol1 = conv_vmf_sol1.mu1 .* mu_r1 + ...
       conv_vmf_sol1.mu2 .* mu_r2 + conv_vmf_sol1.c;
   

mu1 = -conv_const .* rotatedMu_1;
mu2 = -conv_const .* rotatedMu_2;

conv_vmf_sol2.mu1 = R_11 .* mu1 + R_12 .* mu2;
conv_vmf_sol2.mu2 = R_21 .* mu1 + R_22 .* mu2;

% normlize to ideal result
% take absolute, maximal direction is mu_r
[w_max_1,w_max_2] = movmfAbsMu2D(conv_vmf_sol2);

% maximal value in maximal direction
log_estimated_conv_max = conv_vmf_sol2.mu1 .* w_max_1 + ...
                         conv_vmf_sol2.mu2 .* w_max_2;
                     
C = sqrt((apertureVmfL.mu1 + kappaG .* w_max_1).^2 + ...
         (apertureVmfL.mu2 + kappaG .* w_max_2).^2);

c = apertureVmfL.c + scatteringMovmf.c;
log_exact_conv_max = c+log(besseli(0,C,1)) + abs(real(C)) + log(2*pi);

% fix the maximum where it matters
conv_vmf_sol2.c = - log_estimated_conv_max + log_exact_conv_max;
   
sol2 = conv_vmf_sol2.mu1 .* mu_r1 + ...
       conv_vmf_sol2.mu2 .* mu_r2 + conv_vmf_sol2.c;
   
sol1Idx = real(sol1) > real(sol2);

conv_vmf.mu1 = sol1Idx .* conv_vmf_sol1.mu1 + (~sol1Idx) .* conv_vmf_sol2.mu1;
conv_vmf.mu2 = sol1Idx .* conv_vmf_sol1.mu2 + (~sol1Idx) .* conv_vmf_sol2.mu2;
conv_vmf.alpha = conv_vmf_sol1.alpha;
conv_vmf.c = sol1Idx .* conv_vmf_sol1.c + (~sol1Idx) .* conv_vmf_sol2.c;
conv_vmf.dim = conv_vmf_sol1.dim;

end

