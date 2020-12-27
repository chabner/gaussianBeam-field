clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

run('fitting_config.m');
config = preprocess_near_field(config);

plotFitting(config.movmf.cacheIdx);