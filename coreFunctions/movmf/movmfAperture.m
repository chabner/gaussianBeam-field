function [apertureVmf] = movmfAperture(aperuteStd,focalPoints,focalPointsSign,direction,focalPointDim,directionDim)
% Build vmf aperture
%
% INPUT:
% aperuteStd: std of the aperture.
% focalPoints: all the focal points, size of [3,N]
% focalPointsSign: +/- 1
% direction: aperture direction [3,N]
%
% OUTPUT:
% apertureVmf: vmf of the aperture, size of N.

direction_1 = direction(1,:); direction_1 = direction_1(:); direction_1 = shiftdim(direction_1,-directionDim+1);
direction_2 = direction(2,:); direction_2 = direction_2(:); direction_2 = shiftdim(direction_2,-directionDim+1);
direction_3 = direction(3,:); direction_3 = direction_3(:); direction_3 = shiftdim(direction_3,-directionDim+1);

focalPoint_1 = focalPoints(1,:); focalPoint_1 = focalPoint_1(:); focalPoint_1 = shiftdim(focalPoint_1,-focalPointDim+1);
focalPoint_2 = focalPoints(2,:); focalPoint_2 = focalPoint_2(:); focalPoint_2 = shiftdim(focalPoint_2,-focalPointDim+1);
focalPoint_3 = focalPoints(3,:); focalPoint_3 = focalPoint_3(:); focalPoint_3 = shiftdim(focalPoint_3,-focalPointDim+1);

apertureVmf.mu1 = complex(1/(aperuteStd^2) * direction_1 + focalPointsSign * 2 * pi * 1i * focalPoint_1);
apertureVmf.mu2 = complex(1/(aperuteStd^2) * direction_2 + focalPointsSign * 2 * pi * 1i * focalPoint_2);
apertureVmf.mu3 = complex(1/(aperuteStd^2) * direction_3 + focalPointsSign * 2 * pi * 1i * focalPoint_3);

kappa_r = sqrt(real(apertureVmf.mu1).^2 + real(apertureVmf.mu2).^2 + real(apertureVmf.mu3).^2);

apertureVmf.c = ones(size(kappa_r)) * ...
    ((3/2-1)*log(kappa_r(1)) - (3/2)*log(2*pi) - logbesseli(kappa_r(1)));
apertureVmf.alpha = apertureVmf.c.^0;
apertureVmf.c = complex(apertureVmf.c);

apertureVmf.dim = size(apertureVmf.alpha);
apertureVmf.dim(end+1:4) = 1;
apertureVmf.dirDim = directionDim;
end