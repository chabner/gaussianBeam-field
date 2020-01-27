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

if(vmfApperture.dim ~= 3)
    error('working only in 3D')
end

missX = 0;

%% sample first scatterer
Naperture = vmfApperture.N;

% sample (x,y)
kappaMu_r = real(vmfApperture.mu);
kappa_r = sqrt(sum(kappaMu_r.*(conj(kappaMu_r)),3));
P0 = -imag(vmfApperture.mu)/(2*pi);
mu_r = kappaMu_r ./ kappa_r;

mu_r(~isfinite(mu_r(:,:,1:2))) = 0;
mu_r(~isfinite(mu_r(:,:,3))) = 1;

P0 = permute(P0,[1,3,2]);
mu_r = permute(mu_r,[1,3,2]);

while(1)
    % sample one of the apertures
    apertureNum = randi(Naperture);
    
    kappa_g = smpPreprocess(apertureNum).kappa_g;

    icdf = smpPreprocess(apertureNum).icdf;
    NzSamples=length(icdf);
    
    mu_r_n = mu_r(apertureNum,:);
    P0_n = P0(apertureNum,:);
    
    % sample z
    cdfIdx=ceil(rand*NzSamples);
    Pz = smpPreprocess(apertureNum).grid(icdf(cdfIdx));
    
    sigma_gal = sqrt(((kappa_g*mu_r_n(3))^2 + (2*pi*(Pz-P0_n(3))).^2)./(4*pi^2*kappa_g*mu_r_n(3)));
    sigma_hat = sqrt(sigma_gal.^2 + 2*((P0_n(3) - Pz).^2)/kappa_r(1));
    P0_n_gal = P0_n(1:2) - 2 * (mu_r_n(1:2)/mu_r_n(3))*(P0_n(3)-Pz);
    
    Px = normrnd(P0_n_gal(1), sigma_hat);
    Py = normrnd(P0_n_gal(2), sigma_hat);
    
    if(any([Px;Py] < box_min(1:end-1)) || any([Px;Py] > box_max(1:end-1)))
        missX = missX + 1;
        continue;
    end
    
    break;
end

x = [Px;Py;Pz];
px = 0;

for lightNum = 1:1:Naperture
    prob_z = interp1(smpPreprocess(lightNum).grid, smpPreprocess(lightNum).pdf, Pz, 'linear', 1) / ...
        (smpPreprocess(lightNum).grid(2) - smpPreprocess(lightNum).grid(1));
    
    mu_r_n = mu_r(apertureNum,:);
    P0_n = P0(lightNum,:);
    kappa_g = smpPreprocess(apertureNum).kappa_g;
    
    sigma_gal = sqrt(((kappa_g*mu_r_n(3))^2 + (2*pi*(Pz-P0_n(3))).^2)./(4*pi^2*kappa_g*mu_r_n(3)));
    sigma_hat = sqrt(sigma_gal.^2 + 2*((P0_n(3) - Pz).^2)/kappa_r(1));
    P0_n_gal = P0_n(1:2) - 2 * (mu_r_n(1:2)/mu_r_n(3))*(P0_n(3)-Pz);
    
    prob_x = normpdf(Px,P0_n_gal(1),sigma_hat);
    prob_y = normpdf(Py,P0_n_gal(2),sigma_hat);

    px = px + prob_x * prob_y * prob_z;
end

px = px / Naperture;

end

