clear

configFiles = dir(['expirements',filesep,'*.m']);

for configFilesList = 1:1:numel(configFiles)
    clear config
    run([configFiles(configFilesList).folder,filesep,configFiles(configFilesList).name]);
    
    [config,gpuFunc] = preprocessConfig(config);

    Nv = numel(config.focalPointsV.base);
    Nl = numel(config.focalPointsL.base);
    tic                
    [u_nf,us_nf] = run_rendering(config,gpuFunc);
    toc
    u_nf = reshape(u_nf,Nl,Nv,Nv);
    us_nf = reshape(us_nf,Nl,Nv,Nv);
    u_nf = permute(u_nf,[2,3,1]);
    us_nf = permute(us_nf,[2,3,1]);

    config = preprocessConfig_ff(config);
    l = -1:0.02:1;
    v = -1:0.02:1;
    tic
    [u_ff,us_ff] = farField_refocuse(config,l,v);
    toc

    f_m = figure;
    
    imagesc(abs([ ...
        reshape(u_nf, size(u_nf,1) , []) ; ...
        reshape(u_ff, size(u_ff,1) , []) ...
        ]));
    
    filename = configFiles(configFilesList).name;
    filename = filename(1:end-2);
    savefig(f_m,[filename,'_','multipleScattering.fig'])

    f_s = figure;
    
    imagesc(abs([ ...
        reshape(us_nf, size(u_nf,1) , []) ; ...
        reshape(us_ff, size(u_ff,1) , []) ...
        ]));
    
    savefig(f_s,[filename,'_','singleScattering.fig'])
end