function [apertureVmf] = movmfAperture2D(aperuteStd,focalPoints,focalPointsSign,direction,focalPointDim,directionDim,shiftVec,shiftDim)
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

dim = 2;

direction_1 = direction(1,:); direction_1 = direction_1(:); direction_1 = shiftdim(direction_1,-directionDim+1);
direction_2 = direction(2,:); direction_2 = direction_2(:); direction_2 = shiftdim(direction_2,-directionDim+1);

focalPoint_1 = focalPoints(1,:); focalPoint_1 = focalPoint_1(:); focalPoint_1 = shiftdim(focalPoint_1,-focalPointDim+1);
focalPoint_2 = focalPoints(2,:); focalPoint_2 = focalPoint_2(:); focalPoint_2 = shiftdim(focalPoint_2,-focalPointDim+1);

if(nargin > 6)
    shiftPoint_1 = shiftVec(1,:); shiftPoint_1 = shiftPoint_1(:); shiftPoint_1 = shiftdim(shiftPoint_1,-shiftDim+1);
    shiftPoint_2 = shiftVec(2,:); shiftPoint_2 = shiftPoint_2(:); shiftPoint_2 = shiftdim(shiftPoint_2,-shiftDim+1);

    apertureVmf.mu1 = complex(1/(aperuteStd^2) * direction_1 + focalPointsSign * 2 * pi * 1i * (focalPoint_1 + shiftPoint_1));
    apertureVmf.mu2 = complex(1/(aperuteStd^2) * direction_2 + focalPointsSign * 2 * pi * 1i * (focalPoint_2 + shiftPoint_2));
else
    apertureVmf.mu1 = complex(1/(aperuteStd^2) * direction_1 + focalPointsSign * 2 * pi * 1i * focalPoint_1);
    apertureVmf.mu2 = complex(1/(aperuteStd^2) * direction_2 + focalPointsSign * 2 * pi * 1i * focalPoint_2);
end

kappa_r = sqrt(real(apertureVmf.mu1).^2 + real(apertureVmf.mu2).^2);

apertureVmf.c = ones(size(kappa_r)) * ...
    ( - log(2*pi) - logbesseli(dim,kappa_r(1)));
apertureVmf.alpha = apertureVmf.c.^0;
apertureVmf.c = complex(apertureVmf.c);

apertureVmf.dim = size(apertureVmf.alpha);
apertureVmf.dim(end+1:4) = 1;
apertureVmf.dirDim = directionDim;
apertureVmf.focalDim = focalPointDim;
end