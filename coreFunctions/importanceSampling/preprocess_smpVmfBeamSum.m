function [prepocessRes] = preprocess_smpVmfBeamSum(config,apertureVmf_l)
% preprocess the z cdf for convenient sampling

% get geometry / scattering parameters

k = 2*pi;
kappa_r = 1/config.mask_varL^2;
sigt = 1/config.MFP;
kappa_g = config.movmf.mu3;
box_min = config.box_min;
box_max = config.box_max;

% get apertures focal and directions

[mu_r1, mu_r2, mu_r3] = movmfAbsMu(apertureVmf_l);
P0_x = -imag(apertureVmf_l.mu1)/(2*pi);
P0_y = -imag(apertureVmf_l.mu2)/(2*pi);
P0_z = -imag(apertureVmf_l.mu3)/(2*pi);

% project all box points to apertures main axis

% allCubePoints = ...
%     [box_min(1), box_min(1), box_min(1), box_min(1), box_max(1), box_max(1), box_max(1), box_max(1) ; ...
%      box_min(2), box_min(2), box_max(2), box_max(2), box_min(2), box_min(2), box_max(2), box_max(2) ; ...
%      box_min(3), box_max(3), box_min(3), box_max(3), box_min(3), box_max(3), box_min(3), box_max(3)];
% allCubePoints = allCubePoints.';

middlePoint = (box_min + box_max)/2;

projectionX = (middlePoint(1) - P0_x) .* mu_r1;
projectionY = (middlePoint(2) - P0_y) .* mu_r2;
projectionZ = (middlePoint(3) - P0_z) .* mu_r3;

dz = (box_max(3) - box_min(3)) ./ abs(mu_r3);

minZ = projectionX + projectionY + projectionZ - dz/2;
maxZ = projectionX + projectionY + projectionZ + dz/2;

% integrate from minZ to maxZ
Ei = @(x) -expint(-x);
a = sqrt(kappa_r * (kappa_r + kappa_g) / k^2);

if(all(abs(mu_r3(:)) > 0.99999))
    prepocessRes.tau = 0;
else
    prepocessRes.tau = 0.1;
end

l_att = -(1./a) .* exp(-sigt .* (P0_z - minZ)) .* imag(exp(1i .* sigt .* a) .* ...
    (Ei(sigt .* (-1i * a + P0_z - maxZ)) - Ei(sigt .* (-1i * a + P0_z - minZ))) );

l = prepocessRes.tau * (pi ./ a) + (1 - prepocessRes.tau) * l_att;

% sample z on the main aperture axis
% prepocessRes.smp_logpdf = @(g,n,x) (~(x >= minZ(n) & x <= maxZ(n)))*(-1e20) + (-sigt * (x - minZ(n)) - log(a(g).^2 + (x-P0_z(n)).^2));

prepocessRes.smp_cdf = @(g,n,x)  -(1./(l_att((n - 1) * size(l_att,1) + g)*a(g))) .* exp(-sigt .* (P0_z(n) - minZ(n))) .* imag(exp(1i .* sigt .* a(g)) .* ...
    (Ei(sigt .* (-1i * a(g) + P0_z(n) - x)) - Ei(sigt .* (-1i * a(g) + P0_z(n) - minZ(n)))) );

prepocessRes.smp_cdf_a = @(g,n)  -(1./(l_att((n - 1) * size(l_att,1) + g)*a(g))) .* exp(-sigt .* (P0_z(n) - minZ(n)));
prepocessRes.smp_cdf_b = @(g,n)  exp(1i .* sigt .* a(g));
prepocessRes.smp_cdf_c = @(g,n)  Ei(sigt .* (-1i * a(g) + P0_z(n) - minZ(n)));
prepocessRes.smp_cdf_fx = @(g,n,x)  Ei(sigt .* (-1i * a(g) + P0_z(n) - x));

prepocessRes.z0 = @(n) (minZ(n) + maxZ(n))/2;

% calculate the probability to sample point on box
prepocessRes.project = @(x,y,z) ...
    (x - P0_x) .* mu_r1 +       ...
    (y - P0_y) .* mu_r2 +       ...
    (z - P0_z) .* mu_r3;
    
% project to z before usage
prepocessRes.eval_pdf = @(z) 1 ./ (l .* (a.^2 + (z-P0_z).^2)) .* ...
    (prepocessRes.tau + (1 - prepocessRes.tau) .* (z >= minZ & z <= maxZ) .* exp(-sigt .* (z - minZ)) );

prepocessRes.alpha = config.movmf.alpha / sum(config.movmf.alpha);

prepocessRes.w0 = @(g,n,z) sqrt((kappa_r + kappa_g(g)) / k^2 + (z - P0_z(n)).^2 /kappa_r);
prepocessRes.w0_all = @(z) sqrt((kappa_r + kappa_g) / k^2 + (z - P0_z).^2 /kappa_r);

prepocessRes.cauchyInvCdf = @(g,n,x) a(g).*tan(pi * (x - 0.5)) + P0_z(n);

prepocessRes.eval_pdf_att = @(z) 1 ./ (l_att .* (a.^2 + (z-P0_z).^2)) .* ...
    ((z >= minZ & z <= maxZ) .* exp(-sigt .* (z - minZ)) );

prepocessRes.eval_pdf_cauchy = @(z) (a./pi) ./ ((a.^2 + (z-P0_z).^2));

end
