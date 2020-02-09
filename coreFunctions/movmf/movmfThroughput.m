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


x1 = x(1);
x2 = x(2);
x3 = x(3);

dz = dz(:);
Nz = numel(dz);

x1 = repmat(x1,1,Nz);
x2 = repmat(x2,1,Nz);
x3 = repmat(x3,1,Nz);

tmp_throughputVmf.mu1 = signz * 2*pi*1i*x1;
tmp_throughputVmf.mu2 = signz * 2*pi*1i*x2;
tmp_throughputVmf.mu3 = signz * 2*pi*1i*x3;

tmp_throughputVmf.alpha = tmp_throughputVmf.mu1.^0;
tmp_throughputVmf.c = -(sigt/2)*dz.';
tmp_throughputVmf.dim = size(tmp_throughputVmf.alpha);

throughputVmf = movmfMultiple(apertureVmf,tmp_throughputVmf,true);
end


