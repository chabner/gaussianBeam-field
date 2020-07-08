clear
addpath(genpath('.'));
rng('shuffle')

configFiles = dir(['expirements',filesep,'*.m']);

for configFilesList = 1:1:numel(configFiles)
    clear config
    run([configFiles(configFilesList).folder,filesep,configFiles(configFilesList).name]);    
    config = preprocessConfig(config);

    tic
    [C,Cs] = run_rendering(config);

    t_nf = toc
    
    figure, imagesc(abs(reshape(C,size(C,1),[]))), title('NF'), colorbar

    clear config;
    run([configFiles(configFilesList).folder,filesep,configFiles(configFilesList).name]);
    
    config.iterationsRender = config.iterationsRender * config.multiplePaths;
    config.multiplePaths = 1;
    
    config = preprocessConfig_ff(config);

    tic
    [C_ff,Cs_ff] = run_farField(config);
    t_ff = toc
    
    figure, imagesc(abs(reshape(C_ff,size(C_ff,1),[]))), title('FF + refocus'), colorbar
    
%     tic
%     [C_binary,Cs_binary] = farField_refocuse_correlation(config,l,v,focalPointsL,focalPointsL_2,focalPointsV,0,...
%         'forward',[0;0;1],[0;0;1],0.375,true);
%     t_binary = toc
%     
%     figure, imagesc(abs(reshape(C_binary,Nv,Nv*Nl))), title('binary + refocus')

end