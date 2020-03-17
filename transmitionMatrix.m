   
allExpirements = dir(['expirements',filesep,'*.m']);

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
        clear config
        
        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);
        [config,gpuFunc] = preprocessConfig(config);
        Nv = numel(config.focalPointsV.base);
        Nl = numel(config.focalPointsL.base);
        
        tic
        [u,us] = run_rendering(config,gpuFunc);
        u = reshape(u,Nv,Nv,Nl,Nl);
        us = reshape(us,Nv,Nv,Nl,Nl);
        t = toc

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'u','us','config','t');
    end
end
