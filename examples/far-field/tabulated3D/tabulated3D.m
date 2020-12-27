clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

%% hg
run('config_file.m');
config = preprocess_far_field(config);

tic
C_hg = run_far_field(config);
toc

clear config

%% hg tabulated
run('config_file.m');

% in 3d the range of thera is from 0 to pi. in 2d from 0 to 2*pi
theta = linspace(0,pi,1e4);
hg_vals = evaluateHG(theta, config.scatter.g, 3);

config.scatter.type = 2;
config.scatter.f = sqrt(hg_vals);

config = preprocess_far_field(config);

tic
C_tabulated = run_far_field(config);
toc

%% plot

maxVal = max(abs([C_hg(:);C_tabulated(:)]));

figure
subplot(2,4,1)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_hg(:,:,1)))
title(['dl = ',num2str(config.ff.parameters.dl(1))])
ylabel('hg')

subplot(2,4,2)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_hg(:,:,2)))
title(['dl = ',num2str(config.ff.parameters.dl(2))])

subplot(2,4,3)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_hg(:,:,3)))
title(['dl = ',num2str(config.ff.parameters.dl(3))])

subplot(2,4,4)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_hg(:,:,4)))
title(['dl = ',num2str(config.ff.parameters.dl(4))])

subplot(2,4,5)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_tabulated(:,:,1)))
ylabel('tabulated')

subplot(2,4,6)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_tabulated(:,:,2)))

subplot(2,4,7)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_tabulated(:,:,3)))

subplot(2,4,8)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_tabulated(:,:,4)))
