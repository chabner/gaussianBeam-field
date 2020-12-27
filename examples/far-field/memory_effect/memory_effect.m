clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

run('config_file.m');
config = preprocess_far_field(config);
tic
C = run_far_field(config);
toc

plot(config.ff.parameters.dl(:), squeeze(abs(sum(sum(C)))) ./ abs(sum(sum(C(:,:,1)))) )
xlabel('\Delta')
ylabel('C(\Delta)')
title('The memory effect')
