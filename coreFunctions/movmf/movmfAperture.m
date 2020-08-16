function [apertureVmf] = movmfAperture(wavelength, aperuteStd,focalPointsSign,...
    focalPoints_x,focalPoints_y,focalPoints_z,...
    focalDirections_x,focalDirections_y,focalDirections_z)
% Build vmf aperture

dim = 3;

apertureVmf.mu1 = complex(1/(aperuteStd^2) * focalDirections_x + focalPointsSign * 2 * pi * 1i * focalPoints_x./wavelength);
apertureVmf.mu2 = complex(1/(aperuteStd^2) * focalDirections_y + focalPointsSign * 2 * pi * 1i * focalPoints_y./wavelength);
apertureVmf.mu3 = complex(1/(aperuteStd^2) * focalDirections_z + focalPointsSign * 2 * pi * 1i * focalPoints_z./wavelength);

kappa_r = sqrt(real(apertureVmf.mu1).^2 + real(apertureVmf.mu2).^2 + real(apertureVmf.mu3).^2);

apertureVmf.c = ones(size(kappa_r)) * ...
    ((3/2-1)*log(kappa_r(1)) - (3/2)*log(2*pi) - logbesseli(dim,kappa_r(1)));
apertureVmf.c = complex(apertureVmf.c);

apertureVmf.dir1 = focalDirections_x;
apertureVmf.dir2 = focalDirections_y;
apertureVmf.dir3 = focalDirections_z;

end