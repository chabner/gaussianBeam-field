clear
rng('shuffle')
    
%% no importance sampling, g = 0.5
run('configFile');
numOfIterations = config.iterationsRender * config.multiplePaths;

config.g = 0.5;
config.sampleFlag = 11;
config = preprocessConfig(config);

tic
C = run_rendering(config);
t = toc

figure
imagesc(reshape(abs(C),size(C,1),[]))
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
imagesc(reshape(abs(C),size(C,1),[]))
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
imagesc(reshape(abs(C),size(C,1),[]))
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
imagesc(reshape(abs(C),size(C,1),[]))
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
imagesc(reshape(abs(C),size(C,1),[]))
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
imagesc(reshape(abs(C),size(C,1),[]))
colorbar

title('IS, g = 0.9')

clear config