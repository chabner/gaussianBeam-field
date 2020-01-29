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

tmp_throughputVmf.dim = apertureVmf.dim;
tmp_throughputVmf.N = 1;
tmp_throughputVmf.k = 1;

x = x(:);
dz = dz(:);
Nz = numel(dz);
x = repmat(x,1,Nz);

tmp_throughputVmf.mu = permute(signz * 2*pi*1i*(x.'),[1,3,2]);
tmp_throughputVmf.alpha = ones(Nz,1);
tmp_throughputVmf.c = -(sigt/2)*dz;

throughputVmf = movmfMultiple(apertureVmf,tmp_throughputVmf,Nz == 1,true);
end


