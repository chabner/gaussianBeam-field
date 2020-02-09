function [prepocessRes] = preprocess_smpVmfBeamSum(config,zSamplesNum)
% preprocess the z cdf for convenient sampling

sampleGridZ = linspace(config.box_min(3),config.box_max(3),zSamplesNum);

kappa_r = 1/config.mask_varL^2;
sigt = 1/config.attMFP;

kappa_g = config.kappaG;

fp = config.focalPointsL;

prepocessRes.grid = sampleGridZ;
prepocessRes.pdf = zeros(size(fp,2),zSamplesNum);
prepocessRes.icdf = zeros(size(fp,2),zSamplesNum);
prepocessRes.kappa_g = kappa_g;


for lightComp = 1:1:size(fp,2)
    P0_n = fp(:,lightComp);
    P3 = P0_n(3);
    
    mu_n = config.focalDirectionsL(:,lightComp);
    mu3 = mu_n(3);
    
    sigma_gal = sqrt(((kappa_g*mu3)^2 + (2*pi*(sampleGridZ-P3)).^2)./(4*pi^2*kappa_g*mu3));
    sigma_hat = sqrt(sigma_gal.^2 + ((P3 - sampleGridZ).^2)/kappa_r);
    
    Pz = [P0_n(1:2) + 0 * sampleGridZ(:).' ; sampleGridZ(:).'];
    
    dz = cubeDist(Pz,config.box_min,config.box_max,-1 * mu_n);
    
    sampleZfunction = exp(-sigt * dz) .* (sigma_gal./sigma_hat).^2 .* ...
        exp(((kappa_g^2)/(4*pi^2) + ((sampleGridZ - P3)/mu3).^2) .* (1 - mu3^2) ./ sigma_gal.^2);
    
    preprocessPdf = sampleZfunction / sum(sampleZfunction);
    zCdf = cumsum(preprocessPdf);

    preprocessZiCdf = cumsum(hist(zCdf,numel(zCdf)));
    
    prepocessRes.pdf(lightComp,:) = preprocessPdf;
    prepocessRes.icdf(lightComp,:) = preprocessZiCdf;
end

end