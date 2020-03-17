clear

load('is.mat')
expirementType = expirementsResRandPixel1e5;

figure
subplot(5,2,1)
plotExpirement(expirementType.configFile11_g09);
title('g09, random x random w')
colorbar

subplot(5,2,3)
plotExpirement(expirementType.configFile12_g09);
title('g09, exponent x random w')
colorbar

subplot(5,2,5)
plotExpirement(expirementType.configFile14_g09);
title('g09, gaussian sum x random w')
colorbar

subplot(5,2,7)
plotExpirement(expirementType.configFile44_g09);
title('g09, gaussian sum x gaussian sum w, other beam')
colorbar

subplot(5,2,9)
plotExpirement(expirementType.configFile54_g09);
title('g09, gaussian sum x gaussian sum w, same beam')
colorbar

subplot(5,2,2)
plotExpirement(expirementType.configFile11_g05);
title('g05, random x random w')
colorbar

subplot(5,2,4)
plotExpirement(expirementType.configFile14_g05);
title('g05, exponent x random w')
colorbar

subplot(5,2,6)
plotExpirement(expirementType.configFile14_g05);
title('g05, gaussian sum x random w')
colorbar

subplot(5,2,8)
plotExpirement(expirementType.configFile44_g05);
title('g05, gaussian sum x gaussian sum w, other beam')
colorbar

subplot(5,2,10)
plotExpirement(expirementType.configFile54_g05);
title('g05, gaussian sum x gaussian sum w, same beam')
colorbar


function plotExpirement(expStr)
    imagesc(reshape(abs(expStr.C ./ (expStr.itersNum + expStr.xRep)) , size(expStr.C,1) , []) );
end