clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));

run('config_file.m');
config = preprocess_near_field(config);
tic
u = run_near_field(config);
toc

figure
Nl = size(u,3);

for lNum = 1:1:Nl
    subplot(2,ceil(Nl/2),lNum);
    imagesc(config.nf.parameters.v_x,config.nf.parameters.v_y,abs(u(:,:,lNum)))
end
