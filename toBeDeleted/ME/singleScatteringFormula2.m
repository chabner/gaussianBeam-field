clear

taho_grid = -3000:50:3000;
taho_grid = taho_grid*1;
[tahox,tahoy] = ndgrid(taho_grid,taho_grid);

angleNum = 11;

%%
run('config_OD_1_ff')

delta = 1e3;

ox = 0; oy = 0; oz = 0;

mu3 = (1/config.mask_varL)^2;
z0 = 2e5;

%%
alpha1 = sqrt(-4*pi^2 .* ((ox-delta/2).^2 + oy.^2 + (oz - z0).^2) + mu3.^2 - ...
    1i * 4*pi* mu3 * (oz - z0) );

alpha2 = sqrt(-4*pi^2 .* ((ox+delta/2).^2 + oy.^2 + (oz - z0).^2) + mu3.^2 - ...
    1i * 4*pi* mu3 * (oz - z0) );

beta1 = sqrt(-4*pi^2 .* ((ox-delta/2-tahox).^2 + (oy-tahoy).^2 + (oz - z0).^2) + mu3.^2 + ...
    1i * 4*pi* mu3 * (oz - z0) );

beta2 = sqrt(-4*pi^2 .* ((ox+delta/2-tahox).^2 + (oy-tahoy).^2 + (oz - z0).^2) + mu3.^2 + ...
    1i * 4*pi* mu3 * (oz - z0) );

%%
u1 = exp(alpha1 + beta1)./(alpha1.*beta1);
u2 = conj(exp(alpha2 + beta2)./(alpha2.*beta2));

f = figure;
subplot(2,3,1)
imagesc(taho_grid,taho_grid,abs(u1))
title('abs')
colorbar

subplot(2,3,2)
imagesc(taho_grid,taho_grid,real(u1))
title('real')
colorbar

subplot(2,3,3)
imagesc(taho_grid,taho_grid,imag(u1))
title('imag')
colorbar

subplot(2,3,4)
imagesc(taho_grid,taho_grid,abs(u2))
title('abs')
colorbar

subplot(2,3,5)
imagesc(taho_grid,taho_grid,real(u2))
title('real')
colorbar

subplot(2,3,6)
imagesc(taho_grid,taho_grid,imag(u2))
title('imag')
colorbar
