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
Naperture = vmfApperture.dim(2);

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

    icdf = smpPreprocess.icdf(apertureNum,:);
    NzSamples=length(icdf);
    
    mu_r_n1 = mu_r1(apertureNum); mu_r_n2 = mu_r2(apertureNum); mu_r_n3 = mu_r3(apertureNum);
    P0_n_x = P0_x(apertureNum); P0_n_y = P0_y(apertureNum); P0_n_z = P0_z(apertureNum);
    
    % sample z
    cdfIdx=ceil(rand*NzSamples);
    Pz = smpPreprocess.grid(icdf(cdfIdx));
    
    sigma_gal = sqrt(((kappa_g*mu_r_n3)^2 + (2*pi*(Pz-P0_n_z)).^2)./(4*pi^2*kappa_g*mu_r_n3));
    sigma_hat = sqrt(sigma_gal.^2 + ((P0_n_z - Pz).^2)/kappa_r);
    P0_n_gal_x = P0_n_x - 2 * (mu_r_n1/mu_r_n3)*(P0_n_z-Pz);
    P0_n_gal_y = P0_n_y - 2 * (mu_r_n2/mu_r_n3)*(P0_n_z-Pz);
    
    Px = normrnd(P0_n_gal_x, sigma_hat / sqrt(2));
    Py = normrnd(P0_n_gal_y, sigma_hat / sqrt(2));
    
    if(any([Px;Py] < box_min(1:end-1)) || any([Px;Py] > box_max(1:end-1)))
        missX = missX + 1;
        continue;
    end
    
    break;
end

x = [Px;Py;Pz];

prob_z = interp1(smpPreprocess.grid, smpPreprocess.pdf', Pz, 'linear', 1) / ...
    (smpPreprocess.grid(2) - smpPreprocess.grid(1));


sigma_gal = sqrt(((kappa_g*mu_r3).^2 + (2*pi*(Pz-P0_z)).^2)./(4*pi^2*kappa_g*mu_r3));
sigma_hat = sqrt(sigma_gal.^2 + ((P0_z - Pz).^2)/kappa_r);
P0_gal_x = P0_x - 2 * (mu_r1./mu_r3).*(P0_z-Pz);
P0_gal_y = P0_y - 2 * (mu_r2./mu_r3).*(P0_z-Pz);

prob_x = normpdf(Px,P0_gal_x,(sigma_hat/sqrt(2)));
prob_y = normpdf(Py,P0_gal_y,(sigma_hat/sqrt(2)));

px = mean(prob_x .* prob_y .* prob_z);


end

