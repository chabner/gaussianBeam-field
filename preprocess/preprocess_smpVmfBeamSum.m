function [prepocessRes] = preprocess_smpVmfBeamSum(config,zSamplesNum)
% preprocess the z cdf for convenient sampling

apertureVmf_l = config.apertureVmf_l;

sampleGridZ = linspace(config.box_min(3),config.box_max(3),zSamplesNum);
sampleGridZ = sampleGridZ(:);

kappa_r = 1/config.mask_varL^2;
sigt = 1/config.attMFP;

kappa_g = config.kappaG;

prepocessRes.grid = sampleGridZ;
prepocessRes.pdf = zeros([zSamplesNum,apertureVmf_l.dim(2:end)]);
prepocessRes.icdf = zeros([zSamplesNum,apertureVmf_l.dim(2:end)]);
prepocessRes.kappa_g = kappa_g;

[mu_r1, mu_r2, mu_r3] = movmfAbsMu(apertureVmf_l);
P0_x = -imag(apertureVmf_l.mu1)/(2*pi);
P0_y = -imag(apertureVmf_l.mu2)/(2*pi);
P0_z = -imag(apertureVmf_l.mu3)/(2*pi);

sigma_gal = sqrt(((kappa_g*mu_r3).^2 + (2*pi*(sampleGridZ-P0_z)).^2)./(4*pi^2*kappa_g*mu_r3));
sigma_hat = sqrt(sigma_gal.^2 + ((P0_z - sampleGridZ).^2)/kappa_r);

Pz = [P0_x(:)' ; P0_y(:)';P0_z(:)'];
mu_r = [mu_r1(:)' ; mu_r2(:)';mu_r3(:)'];
dz = cubeDist(Pz,config.box_min,config.box_max,-1 * mu_r);
dz = reshape(dz,size(mu_r1));

sampleZfunction = exp(-sigt * dz) .* (sigma_gal./sigma_hat).^2 .* ...
        exp(((kappa_g^2)/(4*pi^2) + ((sampleGridZ - P0_z)./mu_r3).^2) .* (1 - mu_r3.^2) ./ sigma_gal.^2);
    
preprocessPdf = sampleZfunction ./ sum(sampleZfunction,1);
zCdf = cumsum(preprocessPdf,1);

zCdf = reshape(zCdf,zSamplesNum,[]);
preprocessZiCdf = cumsum(hist(zCdf,zSamplesNum),1);
preprocessZiCdf = reshape(preprocessZiCdf,size(preprocessPdf));
    
prepocessRes.pdf = preprocessPdf;
prepocessRes.icdf = preprocessZiCdf;

prepocessRes.icdf(prepocessRes.icdf == 0) = 1;

end