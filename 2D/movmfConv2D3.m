function [conv_vmf]=movmfConv2D3(apertureVmfL,scatteringMovmf,dirv,dirDim)
% convolve aperture with scattering function

% INPUT
% apertureVmf: complex vmf aperture.
% scatteringMovmf: vmf mixture of the scattering function, mu only in
% directions of [0,0,1] or [0,0,-1]
%
% OUTPUT
% conv_vmf: convolution result.
perm = 1:1:5; perm(2) = dirDim; perm(dirDim) = 2;

mu_v1 = permute(dirv(1,:),perm); 
mu_v2 = permute(dirv(2,:),perm); 

kappaG = real(scatteringMovmf.mu2);

s = sqrt((apertureVmfL.mu1 + kappaG .* mu_v1).^2 + (apertureVmfL.mu2 + kappaG .* mu_v2).^2);
conv_const = kappaG ./ s;

conv_vmf.mu1 = complex(conv_const .* apertureVmfL.mu1);
conv_vmf.mu2 = complex(conv_const .* apertureVmfL.mu2);

% normlize to ideal result
% take absolute, maximal direction is mu_r
[w_max_1,w_max_2] = movmfAbsMu2D(conv_vmf);

% maximal value in maximal direction
log_estimated_conv_max = conv_vmf.mu1 .* w_max_1 + ...
                         conv_vmf.mu2 .* w_max_2;
                     
C = sqrt((apertureVmfL.mu1 + kappaG .* w_max_1).^2 + ...
         (apertureVmfL.mu2 + kappaG .* w_max_2).^2);

c = apertureVmfL.c + scatteringMovmf.c;
log_exact_conv_max = c+log(besseli(0,C,1)) + abs(real(C)) + log(2*pi);

% fix the maximum where it matters
conv_vmf.c = - log_estimated_conv_max + log_exact_conv_max;

conv_vmf.dim = size(conv_vmf.mu2);
conv_vmf.alpha = repmat(scatteringMovmf.alpha,[1,conv_vmf.dim(2:end)]);

conv_vmf.dim(end+1:4) = 1;
conv_vmf.focalDim = apertureVmfL.focalDim;

end

