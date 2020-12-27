clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

run('config_file.m');
config = preprocess_near_field(config);
tic
C = run_near_field(config);
toc

plot(config.nf.parameters.delta(:), squeeze(abs(sum(sum(C)))) ./ abs(sum(sum(C(:,:,1)))) )
xlabel('\Delta')
ylabel('C(\Delta)')
title('The memory effect')
