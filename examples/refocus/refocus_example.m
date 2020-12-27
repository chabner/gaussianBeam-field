clear
addpath(genpath(['..',filesep,'..',filesep]));

run('refocus_config.m');
config = preprocess_refocus(config);
tic
C_refocus = run_refocus(config);
toc

clear config
run('refocus_config.m');
config = preprocess_near_field(config);
tic
C_nf = run_near_field(config);
toc

figure
subplot(2,2,1)
imagesc(config.nf.parameters.v_x(:),config.nf.parameters.v_y(:),abs(C_refocus(:,:,1)))
colorbar;
title(['refocus , \Delta = ',num2str(config.nf.parameters.delta(1))])

subplot(2,2,2)
imagesc(config.nf.parameters.v_x(:),config.nf.parameters.v_y(:),abs(C_refocus(:,:,2)))
colorbar;
title(['refocus , \Delta = ',num2str(config.nf.parameters.delta(2))])

subplot(2,2,3)
imagesc(config.nf.parameters.v_x(:),config.nf.parameters.v_y(:),abs(C_nf(:,:,1)))
colorbar;
title(['nf , \Delta = ',num2str(config.nf.parameters.delta(1))])

subplot(2,2,4)
imagesc(config.nf.parameters.v_x(:),config.nf.parameters.v_y(:),abs(C_nf(:,:,2)))
colorbar;
title(['nf , \Delta = ',num2str(config.nf.parameters.delta(2))])
