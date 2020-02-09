clear

load('res_tmp.mat')
expirementType = expirementsRes;

figure
subplot(4,1,1)
plotExpirement(expirementType.configFile11);
colorbar

subplot(4,1,2)
plotExpirement(expirementType.configFile14);
colorbar

subplot(4,1,3)
plotExpirement(expirementType.configFile44);
colorbar

subplot(4,1,4)
plotExpirement(expirementType.configFile44);
colorbar


function plotExpirement(expStr)
    imagesc(reshape(abs(expStr.C ./ (expStr.itersNum + expStr.xRep)) , size(expStr.C,1) , []) );
end