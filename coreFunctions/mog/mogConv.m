function [conv_mog]=mogConv(apertureMog,scatteringMog)
% convolve aperture with scattering function

% INPUT
% apertureVmf: complex vmf aperture.
% scatteringMovmf: vmf mixture of the scattering function, mu only in
% directions of [0,0,1] or [0,0,-1]
%
% OUTPUT
% conv_vmf: convolution result.


sumA1A2 = apertureMog.a + scatteringMog.a;
invSumA1A2 = 1./sumA1A2;

Bdiff1 = apertureMog.b1 - scatteringMog.b1;
Bdiff2 = apertureMog.b2 - scatteringMog.b2;

conv_mog.a = apertureMog.a .* scatteringMog.a .* invSumA1A2;
conv_mog.b1 = scatteringMog.b1 .* invSumA1A2  .* apertureMog.a + ...
              apertureMog.b1 .* invSumA1A2  .* scatteringMog.a;
conv_mog.b2 = scatteringMog.b2 .* invSumA1A2  .* apertureMog.a + ...
              apertureMog.b2 .* invSumA1A2  .* scatteringMog.a;
conv_mog.c = apertureMog.c + scatteringMog.c + ...
    (Bdiff1.^2 + Bdiff2.^2) .* invSumA1A2 / 2 + ...
    log((2*pi) .* invSumA1A2);
conv_mog.alpha = apertureMog.alpha .* scatteringMog.alpha;

conv_mog.dim = size(conv_mog.alpha);
conv_mog.dim(end+1:4) = 1;

end
