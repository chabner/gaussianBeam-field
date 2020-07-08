function [conv_vmf]=movmfConv(apertureVmfL,apertureVmfV,scatteringMovmf)
% convolve aperture with scattering function

% INPUT
% apertureVmf: complex vmf aperture.
% scatteringMovmf: vmf mixture of the scattering function, mu only in
% directions of [0,0,1] or [0,0,-1]
%
% OUTPUT
% conv_vmf: convolution result.
% perm = 1:1:5; perm(2) = dirDim; perm(dirDim) = 2;
% 
% mu_v1 = permute(dirv(1,:,:,:,:),perm); 
% mu_v2 = permute(dirv(2,:,:,:,:),perm); 
% mu_v3 = permute(dirv(3,:,:,:,:),perm);

mu_v1 = apertureVmfV.dir1;
mu_v2 = apertureVmfV.dir2;
mu_v3 = apertureVmfV.dir3;

kappaG = real(scatteringMovmf.mu3);

s = sqrt((apertureVmfL.mu1 + kappaG .* mu_v1).^2 + (apertureVmfL.mu2 + kappaG .* mu_v2).^2 + (apertureVmfL.mu3 + kappaG .* mu_v3).^2);
conv_const = kappaG ./ s;

conv_vmf.mu1 = complex(conv_const .* apertureVmfL.mu1);
conv_vmf.mu2 = complex(conv_const .* apertureVmfL.mu2);
conv_vmf.mu3 = complex(conv_const .* apertureVmfL.mu3);

% % fix the maximum where it matters
% conv_vmf.c = s - conv_const .* (apertureVmfL.mu1 .* mu_v1 + apertureVmfL.mu2 .* mu_v2 + apertureVmfL.mu3 .* mu_v3) + ...
%     log(2*pi ./ s) + apertureVmfL.c + scatteringMovmf.c;

% normlize to ideal result
% take absolute, maximal direction is mu_r
[w_max_1,w_max_2,w_max_3] = movmfAbsMu(conv_vmf);

% maximal value in maximal direction
log_estimated_conv_max = conv_vmf.mu1 .* w_max_1 + ...
                         conv_vmf.mu2 .* w_max_2 + ...
                         conv_vmf.mu3 .* w_max_3;
                     
C = sqrt((apertureVmfL.mu1 + kappaG .* w_max_1).^2 + ...
         (apertureVmfL.mu2 + kappaG .* w_max_2).^2 + ...
         (apertureVmfL.mu3 + kappaG .* w_max_3).^2);

c = apertureVmfL.c + scatteringMovmf.c;
log_exact_conv_max = c+C + log(2*pi) - log(C);

% fix the maximum where it matters
conv_vmf.c = - log_estimated_conv_max + log_exact_conv_max;

end
