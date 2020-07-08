clear

run('isConfig');
config = preprocessConfig(config);
smpNum = 1e4;

tic
[x,~,~,n] = smpVmfBeamSum(config.apertureVmf_l,config.smpPreprocess,config.box_min,config.box_max,smpNum);
toc

figure
scatter(x(1,n==1),x(3,n==1))
hold on

for lNum = 2:1:max(n)
    scatter(x(1,n==lNum),x(3,n==lNum))   
end

