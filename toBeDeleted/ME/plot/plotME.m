% load('res_tmp');

figure
subplot(3,1,1)
hold on
plotMeGaraph(expirementsRes.config_OD_10_nf)
plotMeGaraph(expirementsRes.config_OD_10_middle)
plotMeGaraph(expirementsRes.config_OD_10_ff)
xlim([0,80]);

subplot(3,1,2)
hold on
plotMeGaraph(expirementsRes.config_OD_5_nf)
plotMeGaraph(expirementsRes.config_OD_5_middle)
plotMeGaraph(expirementsRes.config_OD_5_ff)
xlim([0,80]);

subplot(3,1,3)
hold on
plotMeGaraph(expirementsRes.config_OD_1_nf)
plotMeGaraph(expirementsRes.config_OD_1_middle)
plotMeGaraph(expirementsRes.config_OD_1_ff)
xlim([0,80]);

function plotMeGaraph(expirement)
    C = expirement.C;
    x = expirement.config.focalPointsL.base;
    Cres = (squeeze(sum(sum(abs(C),1),2))).^2;
    plot(x,Cres / Cres(1));    
end