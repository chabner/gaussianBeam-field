function [x,px,missX,apertureNum] = smpVmfBeamSum(vmfApperture,smpPreprocess,box_min,box_max)
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
[mu_r1, mu_r2, mu_r3, kappa_r] = movmfAbsMu(vmfApperture);
kappa_r = kappa_r(1);
P0_x = -imag(vmfApperture.mu1)/(2*pi);
P0_y = -imag(vmfApperture.mu2)/(2*pi);
P0_z = -imag(vmfApperture.mu3)/(2*pi);

while(1)
    % sample one of the apertures
    apertureNum = randi(Naperture);
    kappa_g = smpPreprocess.kappa_g;
    
    zIdx = smpPreprocess.iz(apertureNum);
    muIdx = smpPreprocess.imu(apertureNum);

    icdf = smpPreprocess.icdf(:,zIdx,muIdx);
    NzSamples=length(icdf);
    
    mu_r_n1 = mu_r1(apertureNum); mu_r_n2 = mu_r2(apertureNum); mu_r_n3 = mu_r3(apertureNum);
    P0_n_x = P0_x(apertureNum); P0_n_y = P0_y(apertureNum); P0_n_z = P0_z(apertureNum);
    
    % sample z
    cdfIdx=ceil(rand*NzSamples);
    Pz = smpPreprocess.grid(icdf(cdfIdx));
    
    sigma_gal = sqrt(((kappa_g*mu_r_n3)^2 + (2*pi*(Pz-P0_n_z)).^2)./(4*pi^2*kappa_g*mu_r_n3));
    sigma_hat = sqrt(sigma_gal.^2 + ((P0_n_z - Pz).^2)/kappa_r);
    
    P0_n_gal_x = P0_n_x - (mu_r_n1/mu_r_n3)*(P0_n_z-Pz);
    P0_n_gal_y = P0_n_y - (mu_r_n2/mu_r_n3)*(P0_n_z-Pz);
    
    xy = normrnd(0, sigma_hat / sqrt(2),1,2);
    r = sqrt(xy(1).^2 + xy(2).^2);
    
    randDirection = randn(1,3);
    pDir = cross(randDirection,[mu_r_n1,mu_r_n2,mu_r_n3]);
    pDir = pDir ./ sqrt(sum(pDir.^2));
    
    Px = P0_n_gal_x + pDir(1) * r;
    Py = P0_n_gal_y + pDir(2) * r;
    Pz = Pz + pDir(3) * r;
    
    if(any([Px;Py;Pz] < box_min) || any([Px;Py;Pz] > box_max))
        missX = missX + 1;
        continue;
    end
    
    break;
end

x = [Px;Py;Pz];

projectionBase = ...
    (x(1) - P0_x) .* mu_r1 + ...
    (x(2) - P0_y) .* mu_r2 + ...
    (x(3) - P0_z) .* mu_r3;

Pz_projection = mu_r3 .* projectionBase + P0_z;
gridIdx = interp1(smpPreprocess.grid, 1:1:numel(smpPreprocess.grid), Pz_projection, 'NEAREST');

pdfIdx = sub2ind(size(smpPreprocess.pdf), gridIdx(:), smpPreprocess.iz, smpPreprocess.imu);
prob_z = smpPreprocess.pdf(pdfIdx) / (smpPreprocess.grid(2) - smpPreprocess.grid(1));

n_x = mu_r1.*projectionBase + P0_x - Px;
n_y = mu_r2.*projectionBase + P0_y - Py;

% orthogonalDistance = sqrt( (mu_r1.*projectionBase + P0_x - Px).^2 + (mu_r2.*projectionBase + P0_y - Py).^2 + (mu_r3.*projectionBase + P0_z - Pz).^2 );

sigma_gal = sqrt(((kappa_g*mu_r3).^2 + (2*pi*(Pz_projection-P0_z)).^2)./(4*pi^2*kappa_g*mu_r3));
sigma_hat = sqrt(sigma_gal.^2 + ((Pz_projection-P0_z).^2)/kappa_r);

prob_x = normpdf(n_x,0,(sigma_hat/sqrt(2)));
prob_y = normpdf(n_y,0,(sigma_hat/sqrt(2)));

px = mean(prob_x(:) .* prob_y(:) .* prob_z(:));

end

