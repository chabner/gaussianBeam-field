function [prepocessRes] = preprocess_smpVmfBeam(config,lightNum,zSamplesNum)
% preprocess the z cdf for convenient sampling

sampleGridZ = linspace(config.box_min(3),config.box_max(3),zSamplesNum);

kappa_r = 1/config.mask_varL^2;
sigt = 1/config.attMFP;

kappa_g = config.movmf1.kappa;

if(isfield(config,'kappaG'))
    kappa_g = config.kappaG;
end

fp = config.focalPointsL{lightNum};
fd = config.focalDirectionsL{lightNum};

P0_1 = fp(:,1);
P0_2 = fp(:,2);

P3 = (P0_1(3)+P0_2(3))/2;

mu_r_1 = fd(:,1);
mu_r_2 = fd(:,2);

sigma_gal = sqrt(((kappa_g*mu_r_1(3))^2 + (2*pi*(sampleGridZ-P3)).^2)./(4*pi^2*kappa_g*mu_r_1(3)));
x0_gal = @(P0,mu_r) P0(1) - 2 * (mu_r(1)/mu_r(3))*(P3-sampleGridZ);
y0_gal = @(P0,mu_r) P0(2) - 2 * (mu_r(2)/mu_r(3))*(P3-sampleGridZ);
sigma_hat = sqrt(sigma_gal.^2 + 2*((P3 - sampleGridZ).^2)/kappa_r);

sampleZfunction = exp(-sigt * (sampleGridZ - config.box_min(3))) .* ...
    ((sigma_gal.^4) ./ ((kappa_g.*mu_r_1(3)).^2 + (2*pi*(P3 - sampleGridZ)).^2)) .* ...
    exp((((kappa_g^2)/(2*pi^2)) + 2 * ((P3-sampleGridZ)/mu_r_1(3)).^2 ) * (mu_r_1(1)^2+mu_r_1(2)^2) ./ (2*sigma_gal.^2)) .* ...
    (1./sigma_hat).^2 .* ...
    exp(((x0_gal(P0_1,mu_r_1) - x0_gal(P0_2,mu_r_2)).^2 + (y0_gal(P0_1,mu_r_1) - y0_gal(P0_2,mu_r_2)).^2)./(-4*sigma_hat.^2));

preprocessPdf = sampleZfunction / sum(sampleZfunction);
zCdf = cumsum(preprocessPdf);

preprocessZiCdf = cumsum(hist(zCdf,numel(zCdf)));

prepocessRes.grid = sampleGridZ;
prepocessRes.pdf = preprocessPdf;
prepocessRes.icdf = preprocessZiCdf;
prepocessRes.kappa_g = kappa_g;

end