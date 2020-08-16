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

    clear config;
    run([configFiles(configFilesList).folder,filesep,configFiles(configFilesList).name]);
    
    config.iterationsRender = config.iterationsRender * config.multiplePaths;
    config.multiplePaths = 1;
    
    config = preprocessConfig_ff(config);

    tic
    [C_ff,Cs_ff] = run_farField(config);
    t_ff = toc
    
    lambdaTotalNum = size(C,4);
    
    figure, imagesc(abs(reshape(C,size(C,1),[]))), title('NF'), colorbar
    
    for wavlengthNum = 1:1:lambdaTotalNum
        subplot(2,lambdaTotalNum,wavlengthNum)
        imagesc(abs(reshape(C(:,:,:,wavlengthNum),size(C,1),[]))), colorbar
        title(['\lambda = ',num2str(config.parameters.lambda(wavlengthNum))]);
        
        if(wavlengthNum == 1)
            ylabel('NF','FontSize',14)
        end
        
        subplot(2,lambdaTotalNum,wavlengthNum + lambdaTotalNum)
        imagesc(abs(reshape(C_ff(:,:,:,wavlengthNum),size(C_ff,1),[]))), colorbar
        
        if(wavlengthNum == 1)
            ylabel('FF + refocus','FontSize',14)
        end
    end

end