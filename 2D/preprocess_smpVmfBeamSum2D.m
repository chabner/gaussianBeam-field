function [prepocessRes] = preprocess_smpVmfBeamSum2D(config,zSamplesNum,apertureVmf_l)
% preprocess the z cdf for convenient sampling

zDim = 2;
muDim = 3;

sampleGridZ = 2*linspace(config.box_min(end),config.box_max(end),zSamplesNum);
sampleGridZ = sampleGridZ(:);

kappa_r = 1/config.mask_varL^2;
sigt = 1/config.attMFP;

kappa_g = config.kappaG;

prepocessRes.grid = sampleGridZ;
prepocessRes.kappa_g = kappa_g;

[mu_r1, mu_r2] = movmfAbsMu2D(apertureVmf_l);
P0_x = -imag(apertureVmf_l.mu1)/(2*pi);
P0_z = -imag(apertureVmf_l.mu2)/(2*pi);

% prepocessRes.muIdx = size(mu_r1,2) * ((1:size(mu_r1,3))-1) + 1;

[P0_z,~,prepocessRes.iz] = unique(P0_z(:),'stable');
[~,~,prepocessRes.imu] = unique(mu_r1(:) + 1i * mu_r2(:),'stable');
P0_z = shiftdim(P0_z,zDim-1);

P0_x = mean(P0_x(:));
mu_r1 = shiftdim(mu_r1,config.focalDirectionsL.dim-1);
mu_r1 = mu_r1(:,1);
mu_r1 = shiftdim(mu_r1,-muDim+1);

mu_r2 = shiftdim(mu_r2,config.focalDirectionsL.dim-1);
mu_r2 = mu_r2(:,1);
mu_r2 = shiftdim(mu_r2,-muDim+1);

sigma_gal = sqrt(((kappa_g*mu_r2).^2 + (2*pi*(sampleGridZ-P0_z)).^2)./(4*pi^2*kappa_g*abs(mu_r2)));
sigma_hat = sqrt(sigma_gal.^2 + ((P0_z - sampleGridZ).^2)/kappa_r);

Pz = [P0_x(:) * (P0_z(:)').^0 ; P0_z(:)'];
mu_r = [mu_r1(:)' ;mu_r2(:)'];

dimsVec = 1:1:muDim;
dimsVec(2) = muDim;
dimsVec(end) = 2;

mu_r = permute(mu_r,dimsVec);
DZ = permute(sampleGridZ,[2,3,4,1]) - Pz(2,:);
X = Pz(1,:) + DZ .* mu_r1 ;
Z = Pz(2,:) + DZ + 0 * mu_r2;
dz = cubeDist(cat(1,X,Z),config.box_min,config.box_max,-1 * mu_r);
dz = permute(dz,[4,2,3,1]);
dz(dz < 0) = 0;

sampleZfunction = exp(-sigt * dz) .* (sigma_gal./sigma_hat) .* ...
        exp(((kappa_g^2)/(4*pi^2) + ((sampleGridZ - P0_z)./mu_r2).^2) .* (1 - mu_r2.^2) ./ sigma_gal.^2);
    
preprocessPdf = sampleZfunction ./ sum(sampleZfunction,1);
zCdf = cumsum(preprocessPdf,1);

zCdf = reshape(zCdf,zSamplesNum,[]);
preprocessZiCdf = cumsum(hist(zCdf,zSamplesNum));
preprocessZiCdf = reshape(preprocessZiCdf,size(preprocessPdf));
    
prepocessRes.pdf = preprocessPdf;
prepocessRes.icdf = preprocessZiCdf;

prepocessRes.icdf(prepocessRes.icdf == 0) = 1;

end
