function [x,px,missX] = smpVmfBeam(vmfApperture,smpPreprocess,box_min,box_max)
% Sample first scatterer and first scattering direction

% INPUT:
% vmfApperture - vmf apperture, must recive 1 or 2 appertures.
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

if(vmfApperture.N ~= 2)
    error('working only for 2 apertures')
end

missX = 0;

%% sample first scatterer

% sample (x,y)
kappaMu = vmfApperture.kappa .* vmfApperture.mu;
kappaMu_r = real(kappaMu);
kappa_r = sqrt(sum(kappaMu_r.*(conj(kappaMu_r)),3));
P0 = -imag(kappaMu)/(2*pi);
mu_r = kappaMu_r ./ kappa_r;

mu_r = mu_r(1,1,:);
P0_1 = P0(1,:);
P0_2 = P0(2,:);

kappa_g = smpPreprocess.kappa_g;

icdf = smpPreprocess.icdf;
N=length(icdf);

while(1)
    % sample z
    cdfIdx=ceil(rand*N);
    Pz = smpPreprocess.grid(icdf(cdfIdx));
    
    sigma_gal = sqrt(((kappa_g*mu_r(3))^2 + (2*pi*(Pz-P0_1(3))).^2)./(4*pi^2*kappa_g*mu_r(3)));
    sigma_hat = sqrt(sigma_gal.^2 + 2*((P0_1(3) - Pz).^2)/kappa_r(1));
    P0_1_gal = P0_1(1:2) - 2 * (mu_r(1:2)/mu_r(3))*(P0_1(3)-Pz);
    P0_2_gal = P0_2(1:2) - 2 * (mu_r(1:2)/mu_r(3))*(P0_2(3)-Pz);
   
    Px = normrnd((P0_1_gal(1) + P0_2_gal(1))/2, sigma_hat);
    Py = normrnd((P0_1_gal(2) + P0_2_gal(2))/2, sigma_hat);
    
    if(any([Px;Py] < box_min(1:end-1)) || any([Px;Py] > box_max(1:end-1)))
        missX = missX + 1;
        continue;
    end
    
    break;
end

x = [Px;Py;Pz];

prob_z = interp1(smpPreprocess.grid, smpPreprocess.pdf, Pz, 'linear', 1) / ...
    (smpPreprocess.grid(2) - smpPreprocess.grid(1));
prob_x = normpdf(Px,(P0_1_gal(1) + P0_2_gal(1))/2,sigma_hat);
prob_y = normpdf(Py,(P0_1_gal(1) + P0_2_gal(1))/2,sigma_hat);

% xProb = normcdf(box_max(1),(P0_1_gal(1) + P0_2_gal(1))/2,sigma_hat) - ...
%     normcdf(box_min(1),(P0_1_gal(1) + P0_2_gal(1))/2,sigma_hat);
% 
% yProb = normcdf(box_max(2),(P0_1_gal(2) + P0_2_gal(2))/2,sigma_hat) - ...
%     normcdf(box_min(2),(P0_1_gal(2) + P0_2_gal(2))/2,sigma_hat);

xProb = 1;
yProb = 1;

px = (prob_x * prob_y * prob_z) / (xProb * yProb);

end

