clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));
load('scattering_functions.mat');

%%
run('config_file.m');
config.simulation.cbs = false;
config.scatter.f = scattering_function_1;
config = preprocess_far_field(config);

tic
C_sct1_nocbs = run_far_field(config);
toc

clear config

%%
run('config_file.m');
config.simulation.cbs = true;
config.scatter.f = scattering_function_1;
config = preprocess_far_field(config);

tic
C_sct1_cbs = run_far_field(config);
toc

clear config

%%
run('config_file.m');
config.simulation.cbs = false;
config.scatter.f = scattering_function_2;
config = preprocess_far_field(config);

tic
C_sct2_nocbs = run_far_field(config);
toc

clear config

%%
run('config_file.m');
config.simulation.cbs = true;
config.scatter.f = scattering_function_2;
config = preprocess_far_field(config);

tic
C_sct2_cbs = run_far_field(config);
toc

clear config

%% plot
run('config_file.m');

f = figure;
f.Position = [0,0,1200,800];
subplot(2,3,1)
plot(view_angles, abs(scattering_function_1).^2,'linewidth',2)
ylabel('phase function 1')
title('phase function')

subplot(2,3,2)
plot(config.ff.parameters.theta_v, abs(C_sct1_nocbs(:,:,1)),'linewidth',2)
title('speckle correlation: no cbs')

subplot(2,3,3)
plot(config.ff.parameters.theta_v, abs(C_sct1_cbs(:,:,1)),'linewidth',2)
title('speckle correlation: cbs')

subplot(2,3,4)
plot(view_angles, abs(scattering_function_2).^2,'linewidth',2)
ylabel('phase function 2')

subplot(2,3,5)
plot(config.ff.parameters.theta_v, abs(C_sct2_nocbs(:,:,1)),'linewidth',2)

subplot(2,3,6)
plot(config.ff.parameters.theta_v, abs(C_sct2_cbs(:,:,1)),'linewidth',2)

sgtitle('forward configuration')

%% plot

f = figure;
f.Position = [0,0,1200,800];
subplot(2,3,1)
plot(view_angles, abs(scattering_function_1).^2,'linewidth',2)
ylabel('phase function 1')
title('phase function')

subplot(2,3,2)
plot(config.ff.parameters.theta_v, abs(C_sct1_nocbs(:,:,2)),'linewidth',2)
title('speckle correlation: no cbs')

subplot(2,3,3)
plot(config.ff.parameters.theta_v, abs(C_sct1_cbs(:,:,2)),'linewidth',2)
title('speckle correlation: cbs')

subplot(2,3,4)
plot(view_angles, abs(scattering_function_2).^2,'linewidth',2)
ylabel('phase function 2')

subplot(2,3,5)
plot(config.ff.parameters.theta_v, abs(C_sct2_nocbs(:,:,2)),'linewidth',2)

subplot(2,3,6)
plot(config.ff.parameters.theta_v, abs(C_sct2_cbs(:,:,2)),'linewidth',2)

sgtitle('backward configuration')