function [apertureVmf] = movmfAperture(aperuteStd,focalPoints,focalPointsSign,direction)
% Build vmf aperture
%
% INPUT:
% aperuteStd: std of the aperture.
% focalPoints: all the focal points, size of [dim,N]
% focalPointsSign: +/- 1
% direction: aperture direction [dim,N]
%
% OUTPUT:
% apertureVmf: vmf of the aperture, size of N.

N = size(focalPoints,2);
dim = size(focalPoints,1);

if(nargin < 4)
    direction = zeros(dim,N);
    direction(end,:) = 1;
end

apertureVmf.dim = dim;
apertureVmf.N = N;
apertureVmf.k = 1;

mu = 1/(aperuteStd^2) * direction.' + focalPointsSign * 2 * pi * 1i * focalPoints.';
apertureVmf.mu = permute(mu,[1,3,2]);

mu_r = real(mu);
kappa_r = sqrt(sum(abs(mu_r).^2,2));
% apertureVmf.c = (dim/2-1)*log(kappa_r) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,kappa_r);
apertureVmf.c = ones(size(kappa_r)) * ...
    ((dim/2-1)*log(kappa_r(1)) - (dim/2)*log(2*pi) - logbesseli(kappa_r(1)));
apertureVmf.alpha = ones(N,1,class(focalPoints));

end