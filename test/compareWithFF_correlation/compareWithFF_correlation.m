clear
addpath(genpath('.'));
rng('shuffle')

configFiles = dir(['expirements',filesep,'*.m']);

for configFilesList = 1:1:numel(configFiles)
    clear config
    run([configFiles(configFilesList).folder,filesep,configFiles(configFilesList).name]);    
    config = preprocessConfig(config);

    Nv = numel(config.focalPointsV.base);
    Nl = numel(config.focalPointsL.base);

    focalPointsL = config.focalPointsL.vector;
    focalPointsV = config.focalPointsV.vector;
    
    focalPointsL_2 = config.focalPointsL.vector_2;
    focalPointsV_2 = config.focalPointsV.vector_2;

    C_nf = 0;
    Cs_nf = 0;

    tic
    [C,Cs] = run_rendering(config);
    C = reshape(C,Nv,Nv,Nl);
    Cs = reshape(C,Nv,Nv,Nl);

    t_nf = toc
    
    figure, imagesc(abs(reshape(C,Nv,Nv*Nl))), title('NF'), colorbar

    clear config;
    run([configFiles(configFilesList).folder,filesep,configFiles(configFilesList).name]);
    config = preprocessConfig_ff(config);
    l = -1:0.1:1;
    v = -1:0.1:1;

    tic
    [C_ff,Cs_ff] = farField_refocuse_correlation(config,l,v,focalPointsL,focalPointsL_2,focalPointsV,0,...
        'forward',[0;0;1],[0;0;1]);
    t_ff = toc
    
    figure, imagesc(abs(reshape(C_ff,Nv,Nv*Nl))), title('FF + refocus'), colorbar
    
    tic
    [C_binary,Cs_binary] = farField_refocuse_correlation(config,l,v,focalPointsL,focalPointsL_2,focalPointsV,0,...
        'forward',[0;0;1],[0;0;1],0.375,true);
    t_binary = toc
    
    figure, imagesc(abs(reshape(C_binary,Nv,Nv*Nl))), title('binary + refocus')

end