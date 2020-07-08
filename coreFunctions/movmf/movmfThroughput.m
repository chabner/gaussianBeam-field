function [throughputVmf] = movmfThroughput(apertureVmf,x,signz,sigt,dz)
% compute the throughput due to scattering
%
% INPUT:
% apertureVmf: the vmf of the aperture
% x: focal point, [dim,1]
% dz: distance from the box edge
% sigt: 1/MFP
% % focalPointsSign: +/- 1
%
% OUTPUT
% throughputVmf: the throughput after scattering

throughputVmf.mu1 = apertureVmf.mu1 - signz * 2*pi*1i*x(1,1,:);
throughputVmf.mu2 = apertureVmf.mu2 - signz * 2*pi*1i*x(2,1,:);
throughputVmf.mu3 = apertureVmf.mu3 - signz * 2*pi*1i*x(3,1,:);
throughputVmf.c   = apertureVmf.c   + complex(-(sigt/2)*dz);

end


