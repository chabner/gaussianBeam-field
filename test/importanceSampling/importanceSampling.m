clear
rng('shuffle')

run('configFile');
Nl = numel(config.focalPointsL.base);
Nv = numel(config.focalPointsV.base);
numOfIterations = config.iterationsRender * config.multiplePaths;
    
%% no importance sampling, g = 0.5
run('configFile');
config.g = 0.5;
config.sampleFlag = 11;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),[Nv,Nv*Nl]))
colorbar

title('no IS, g = 0.5')

clear config


%% position importance sampling, g = 0.5
run('configFile');
config.g = 0.5;
config.sampleFlag = 14;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),[Nv,Nv*Nl]))
colorbar

title('position IS, g = 0.5')

clear config

%% importance sampling, g = 0.5
run('configFile');
config.g = 0.5;
config.sampleFlag = 54;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),[Nv,Nv*Nl]))
colorbar

title('IS, g = 0.5')

clear config

%% no importance sampling, g = 0.9
run('configFile');
config.g = 0.9;
config.sampleFlag = 11;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),[Nv,Nv*Nl]))
colorbar

title('no IS, g = 0.9')

clear config


%% position importance sampling, g = 0.9
run('configFile');
config.g = 0.9;
config.sampleFlag = 14;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),[Nv,Nv*Nl]))
colorbar

title('position IS, g = 0.9')

clear config

%% importance sampling, g = 0.9
run('configFile');
config.g = 0.9;
config.sampleFlag = 54;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),[Nv,Nv*Nl]))
colorbar

title('IS, g = 0.9')

clear config