clear

load('is.mat')
expirementType = expirementsResRandPixel1e5;

figure
subplot(4,1,1)
plotExpirement(expirementType.configFile11_g09);
colorbar

subplot(4,1,2)
plotExpirement(expirementType.configFile14_g09);
colorbar

subplot(4,1,3)
plotExpirement(expirementType.configFile44_g09);
colorbar

subplot(4,1,4)
plotExpirement(expirementType.configFile54_g09);
colorbar

figure
subplot(4,1,1)
plotExpirement(expirementType.configFile11_g05);
colorbar

subplot(4,1,2)
plotExpirement(expirementType.configFile14_g05);
colorbar

subplot(4,1,3)
plotExpirement(expirementType.configFile44_g05);
colorbar

subplot(4,1,4)
plotExpirement(expirementType.configFile54_g05);
colorbar


function plotExpirement(expStr)
    imagesc(reshape(abs(expStr.C ./ (expStr.itersNum + expStr.xRep)) , size(expStr.C,1) , []) );
end