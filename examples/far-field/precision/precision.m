clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

%% single
run('config_file.m');
config = preprocess_far_field(config);

tic
C_single = run_far_field(config);
toc

clear config

%% double
run('config_file.m');
config.simulation.precision = 'double';

config = preprocess_far_field(config);

tic
C_double = run_far_field(config);
toc

%% plot

figure
subplot(1,2,1)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_single(:,:,2)))
title('single')

subplot(1,2,2)
imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(C_double(:,:,2)))
title('double')
