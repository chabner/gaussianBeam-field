clear

figure
load('res_tmp_1.mat')
expirementType = expirementsRes;
subplot(3,1,1)
imagesc(reshape(abs(expirementType.configFile11.C) ,  ...
    size(expirementType.configFile11.C,1) , []) ./ ...
    (expirementType.configFile11.itersNum + expirementType.configFile11.xRep));

load('res_tmp.mat')
expirementType = expirementsRes;

subplot(3,1,2)
imagesc(reshape(abs(expirementType.configFile44.C) ,  ...
    size(expirementType.configFile44.C,1) , []) ./ ...
    (expirementType.configFile44.itersNum + expirementType.configFile44.xRep));

subplot(3,1,3)
imagesc(reshape(abs(expirementType.configFile14.C) ,  ...
    size(expirementType.configFile14.C,1) , []) ./ ...
    (expirementType.configFile14.itersNum + expirementType.configFile14.xRep));

