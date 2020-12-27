clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

run('config_file.m');
config = preprocess_far_field(config);
tic
u = run_far_field(config);
toc

figure
Nl = size(u,3);

for lNum = 1:1:Nl
    subplot(2,ceil(Nl/2),lNum);
    imagesc(config.ff.parameters.vx,config.ff.parameters.vy,abs(u(:,:,lNum)))
end
