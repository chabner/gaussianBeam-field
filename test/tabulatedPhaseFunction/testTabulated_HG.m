clear
g = 0.95;

%% HG
run('configFile.m');
config.g = g;
config.forwardWeight = 1;
config.sctType = 3;

config = preprocessConfig(config);

plotFitting(config.cacheIdx)

Nl = numel(config.focalPointsL.base);
Nv = numel(config.focalPointsV.base);

tic
[C_hg,Cs_hg,xRep] = run_rendering(config);
t = toc

figure, imagesc(abs(reshape(C_hg,Nv,Nv*Nl))), colorbar

%% Tabulated

clear config
run('configFile.m');

config.sctType = 2;

ampfunc.g = g;
ampfunc.forwardWeight = 1;
theta = linspace(0,pi,1e3);

config.pdf = sqrt(evaluateHG(theta, ampfunc, 0, 3));

config = preprocessConfig(config);
plotFitting(config.cacheIdx)

Nl = numel(config.focalPointsL.base);
Nv = numel(config.focalPointsV.base);

tic
[C_tabulated,Cs_tabulated] = run_rendering(config);
t = toc

figure, imagesc(abs(reshape(C_tabulated,Nv,Nv*Nl))), colorbar