function [x,px,missX,apertureNum] = smpVmfBeamSum2D(vmfApperture,smpPreprocess,box_min,box_max)
% Sample first scatterer and first scattering direction

% INPUT:
% vmfApperture - vmf apperture
% sigma_t - 1/MFP
% box_min,box_max - [1,dim], the box size
% movmfScatterer - movmf of the scattering phsae function, size of [1,k]

% OUTPUT:
% x - [1,dim] - first scatterer position
% px - fist scatter position probability
% missX - number of trials failed in x

missX = 0;

%% sample first scatterer
Naperture = prod(vmfApperture.dim(2:end));

% sample (x,y)
[mu_r1, mu_r2, kappa_r] = movmfAbsMu2D(vmfApperture);
kappa_r = kappa_r(1);
P0_x = -imag(vmfApperture.mu1)/(2*pi);
P0_z = -imag(vmfApperture.mu2)/(2*pi);

while(1)
    % sample one of the apertures
    apertureNum = randi(Naperture);
    kappa_g = smpPreprocess.kappa_g;
    
    zIdx = smpPreprocess.iz(apertureNum);
    muIdx = smpPreprocess.imu(apertureNum);

    icdf = smpPreprocess.icdf(:,zIdx,muIdx);
    NzSamples=length(icdf);
    
    mu_r_n1 = mu_r1(apertureNum); mu_r_n2 = mu_r2(apertureNum);
    P0_n_x = P0_x(apertureNum); P0_n_z = P0_z(apertureNum);
    
    % sample z
    cdfIdx=ceil(rand*NzSamples);
    Pz = smpPreprocess.grid(icdf(cdfIdx));
    
    sigma_gal = sqrt(((kappa_g*mu_r_n2)^2 + (2*pi*(Pz-P0_n_z)).^2)./(4*pi^2*kappa_g*abs(mu_r_n2)));
    sigma_hat = sqrt(sigma_gal.^2 + ((P0_n_z - Pz).^2)/kappa_r);
    
    P0_n_gal_x = P0_n_x - (mu_r_n1/mu_r_n2)*(P0_n_z-Pz);
    
    r = normrnd(0, sigma_hat / sqrt(2),1);
    
    randDirection = 2*(rand < 0.5) - 1;
    pDir = randDirection * [mu_r_n2, -mu_r_n1];
    pDir = pDir ./ sqrt(sum(pDir.^2));
    
    Px = P0_n_gal_x + pDir(1) * r;
    Pz = Pz + pDir(2) * r;
    
    if(any([Px;Pz] < box_min) || any([Px;Pz] > box_max))
        missX = missX + 1;
        continue;
    end
    
    break;
end

x = [Px;Pz];

% CHECK!!
projectionBase = ...
    (x(1) - P0_x) .* mu_r1 + ...
    (x(2) - P0_z) .* mu_r2;

Pz_projection = mu_r2 .* projectionBase + P0_z;
gridIdx = interp1(smpPreprocess.grid, 1:1:numel(smpPreprocess.grid), Pz_projection, 'NEAREST');

pdfIdx = sub2ind(size(smpPreprocess.pdf), gridIdx(:), smpPreprocess.iz, smpPreprocess.imu);
prob_z = smpPreprocess.pdf(pdfIdx) / (smpPreprocess.grid(2) - smpPreprocess.grid(1));

n_x = mu_r1.*projectionBase + P0_x - Px;

% orthogonalDistance = sqrt( (mu_r1.*projectionBase + P0_x - Px).^2 + (mu_r2.*projectionBase + P0_y - Py).^2 + (mu_r3.*projectionBase + P0_z - Pz).^2 );

sigma_gal = sqrt(((kappa_g*mu_r2).^2 + (2*pi*(Pz_projection-P0_z)).^2)./(4*pi^2*kappa_g*abs(mu_r2)));
sigma_hat = sqrt(sigma_gal.^2 + ((Pz_projection-P0_z).^2)/kappa_r);

prob_x = normpdf(n_x,0,(sigma_hat/sqrt(2)));

px = mean(prob_x(:) .* prob_z(:));

end

