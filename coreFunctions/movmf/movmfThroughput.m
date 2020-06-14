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

x1 = x(1,1,1,1,1,:);
x2 = x(2,1,1,1,1,:);
x3 = x(3,1,1,1,1,:);

zDim = size(dz);
zDim(end+1:1:6) = 1;
Nz = prod(zDim(1:4));

dimVec = ones(1,6);
dimVec(5:6) = zDim(5:6);
dimVec(apertureVmf.dirDim) = Nz;

dz = reshape(dz,dimVec);
x1 = repmat(x1,[dimVec(1:end-1),1]);
x2 = repmat(x2,[dimVec(1:end-1),1]);
x3 = repmat(x3,[dimVec(1:end-1),1]);

tmp_throughputVmf.mu1 = signz * 2*pi*1i*x1;
tmp_throughputVmf.mu2 = signz * 2*pi*1i*x2;
tmp_throughputVmf.mu3 = signz * 2*pi*1i*x3;

tmp_throughputVmf.alpha = real(tmp_throughputVmf.mu1.^0);
tmp_throughputVmf.c = complex(-(sigt/2)*dz);
tmp_throughputVmf.dim = size(tmp_throughputVmf.alpha);
tmp_throughputVmf.dim(end+1:4) = 1;

if(numel(tmp_throughputVmf.dim) == 6)
    tmp_throughputVmf.dim(6) = [];
end
tmp_throughputVmf.dim(end+1:1:5) = 1;

throughputVmf = movmfMultiple(apertureVmf,tmp_throughputVmf,true);
end


