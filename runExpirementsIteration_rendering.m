    addpath(genpath('../../GBcode'))
    
    allExpirements = dir(['expirements',filesep,'*.m']);
    
    rng('shuffle')
    for iterNum = 1:1:1e3
        for expirementFile = 1:1:numel(allExpirements)
            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

            config = buildConfig_rendering(config);
            tic
            [u_small,u_big,us_small,us_big] = run_rendering(config);
            t = toc

            T = datetime('now','Format','ddMM_HHmmss_SSS');
            save(['res2',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
                'u_small','u_big','config','t','us_small','us_big');
        end
    end
    
    rmpath(genpath('../../GBcode'))
