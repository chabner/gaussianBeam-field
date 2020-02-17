run('isConfig');
config = preprocessConfig(config);
smpNum = 1e5;

x = zeros(3,smpNum);
n = zeros(1,smpNum);
for sNum = 1:1:smpNum
    [x(:,sNum),~,~,n(sNum)] = smpVmfBeamSum(config.apertureVmf_l,config.smpPreprocess,config.box_min,config.box_max);
end

figure
scatter(x(1,n==1),x(3,n==1))
hold on
scatter(x(1,n==2),x(3,n==2))
scatter(x(1,n==3),x(3,n==3))
