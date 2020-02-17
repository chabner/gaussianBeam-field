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

x = complex(x);

x1 = x(1);
x2 = x(2);
x3 = x(3);

Nz = numel(dz);
dimVec = ones(1,apertureVmf.dirDim);
dimVec(end) = Nz;

dz = reshape(dz,dimVec);
x1 = repmat(x1,dimVec);
x2 = repmat(x2,dimVec);
x3 = repmat(x3,dimVec);

tmp_throughputVmf.mu1 = signz * 2*pi*1i*x1;
tmp_throughputVmf.mu2 = signz * 2*pi*1i*x2;
tmp_throughputVmf.mu3 = signz * 2*pi*1i*x3;

tmp_throughputVmf.alpha = real(tmp_throughputVmf.mu1.^0);
tmp_throughputVmf.c = complex(-(sigt/2)*dz);
tmp_throughputVmf.dim = size(tmp_throughputVmf.alpha);
tmp_throughputVmf.dim(end+1:4) = 1;

throughputVmf = movmfMultiple(apertureVmf,tmp_throughputVmf,true);
end


