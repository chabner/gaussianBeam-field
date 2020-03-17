clear

tahox = 0;
tahoy = 500;

angleNum = 11;
x_grid = (-2e4:1e2:2e4)/2;
z_grid = 0;

%%
run('config_OD_1_ff')

Nl = numel(config.focalPointsL.base);
delta = config.focalPointsL.base(angleNum + Nl/2) - config.focalPointsL.base(angleNum);

[ox,oy,oz] = ndgrid(x_grid,x_grid,z_grid);

mu3 = (1/config.mask_varL)^2;
z0 = config.focalPointsL.plain;

%%
alpha1 = sqrt(-4*pi^2 .* ((ox-delta/2).^2 + oy.^2 + (oz - z0).^2) + mu3.^2 - ...
    1i * 4*pi* mu3 * (oz - z0) );

alpha2 = sqrt(-4*pi^2 .* ((ox+delta/2).^2 + oy.^2 + (oz - z0).^2) + mu3.^2 - ...
    1i * 4*pi* mu3 * (oz - z0) );

beta1 = sqrt(-4*pi^2 .* ((ox-delta/2-tahox).^2 + (oy-tahoy).^2 + (oz - z0).^2) + mu3.^2 + ...
    1i * 4*pi* mu3 * (oz - z0) );

beta2 = sqrt(-4*pi^2 .* ((ox+delta/2-tahox).^2 + (oy-tahoy).^2 + (oz - z0).^2) + mu3.^2 + ...
    1i * 4*pi* mu3 * (oz - z0) );

% ZZ = 1i * 2*pi * abs(oz - z0);
% 
% alpha1 = ZZ + 2 * (( -4*pi^2 .* ((ox-delta/2).^2 + oy.^2) ) - 1i * 4*pi* mu3 * (oz - z0)) ./ ZZ - ...
%     0.125 * ((( -4*pi^2 .* ((ox-delta/2).^2 + oy.^2) ) - 1i * 4*pi* mu3 * (oz - z0)) ./ ZZ);
% alpha2 = ZZ + 2 * (( -4*pi^2 .* ((ox+delta/2).^2 + oy.^2) ) - 1i * 4*pi* mu3 * (oz - z0)) ./ ZZ;
% beta1 = ZZ + 2 * (( -4*pi^2 .* ((ox-delta/2-tahox).^2 + (oy-tahoy).^2) ) + 1i * 4*pi* mu3 * (oz - z0)) ./ ZZ;
% beta2 = ZZ + 2 * (( -4*pi^2 .* ((ox+delta/2-tahox).^2 + (oy-tahoy).^2) ) + 1i * 4*pi* mu3 * (oz - z0)) ./ ZZ;


% alpha1 = sqrt(-4*pi^2 .* ((ox-delta/2).^2 + oy.^2 + (oz - z0).^2) - ...
%     1i * 4*pi* mu3 * (oz - z0) );
% 
% alpha2 = sqrt(-4*pi^2 .* ((ox+delta/2).^2 + oy.^2 + (oz - z0).^2) - ...
%     1i * 4*pi* mu3 * (oz - z0) );
% 
% beta1 = sqrt(-4*pi^2 .* ((ox-delta/2-tahox).^2 + (oy-tahoy).^2 + (oz - z0).^2) + ...
%     1i * 4*pi* mu3 * (oz - z0) );
% 
% beta2 = sqrt(-4*pi^2 .* ((ox+delta/2-tahox).^2 + (oy-tahoy).^2 + (oz - z0).^2) + ...
%     1i * 4*pi* mu3 * (oz - z0) );


% alpha1 = sqrt((-1i * 2*pi * (ox - delta/2)).^2 + ...
%     (-1i * 2*pi * oy).^2 + ...
%     (mu3 - 1i * 2*pi * (oz - z0)).^2);
% 
% alpha2 = sqrt((-1i * 2*pi * (ox + delta/2)).^2 + ...
%     (-1i * 2*pi * oy).^2 + ...
%     (mu3 - 1i * 2*pi * (oz - z0)).^2);
% 
% beta1 = sqrt((1i * 2*pi * (ox - delta/2 - tahox)).^2 + ...
%     (1i * 2*pi * (oy - tahoy)).^2 + ...
%     (mu3 + 1i * 2*pi * (oz - z0)).^2);
% 
% beta2 = sqrt((1i * 2*pi * (ox + delta/2 - tahox)).^2 + ...
%     (1i * 2*pi * (oy - tahoy)).^2 + ...
%     (mu3 + 1i * 2*pi * (oz - z0)).^2);

%%
Cs_formula = exp(alpha1 + conj(alpha2) + beta1 + conj(beta2));

f = figure;
subplot(3,1,1)
imagesc(x_grid,x_grid,abs(Cs_formula))
title('abs')
colorbar

subplot(3,1,2)
imagesc(x_grid,x_grid,real(Cs_formula))
title('real')
colorbar

subplot(3,1,3)
imagesc(x_grid,x_grid,imag(Cs_formula))
title('imag')
colorbar

f.Position = [0,0,800,1200];

suptitle(['tahox = ',num2str(round(tahox)),'. tahoy = ',num2str(round(tahoy))])

%% Gaussian

gl1_s = 2*pi ./ (mu3 - 1i * 2*pi * (oz-z0));
gl1_c = ((2*pi*1i*(ox+delta/2)).^2 + (2*pi*1i*oy).^2) ./ (2 * (mu3 - 1i * 2*pi * (oz-z0)));

gl2_s = 2*pi ./ (mu3 - 1i * 2*pi * (oz-z0));
gl2_c = ((2*pi*1i*(ox-delta/2)).^2 + (2*pi*1i*oy).^2) ./ (2 * (mu3 - 1i * 2*pi * (oz-z0)));

gv1_s = 2*pi ./ (mu3 + 1i * 2*pi * (oz-z0));
gv1_c = ((2*pi*1i*(ox+delta/2-tahox)).^2 + (2*pi*1i*(oy-tahoy)).^2) ./ (2 * (mu3 + 1i * 2*pi * (oz-z0)));

gv2_s = 2*pi ./ (mu3 + 1i * 2*pi * (oz-z0));
gv2_c = ((2*pi*1i*(ox-delta/2-tahox)).^2 + (2*pi*1i*(oy-tahoy)).^2) ./ (2 * (mu3 + 1i * 2*pi * (oz-z0)));

Cs_g = exp(gl1_c + conj(gl2_c) + gv1_c + conj(gv2_c) + ...
    log(gl1_s .* conj(gl2_s) .* gv1_s .* conj(gv2_s)));

f = figure;
subplot(3,1,1)
imagesc(x_grid,x_grid,abs(Cs_g))
title('abs')
colorbar

subplot(3,1,2)
imagesc(x_grid,x_grid,real(Cs_g))
title('real')
colorbar

subplot(3,1,3)
imagesc(x_grid,x_grid,imag(Cs_g))
title('imag')
colorbar