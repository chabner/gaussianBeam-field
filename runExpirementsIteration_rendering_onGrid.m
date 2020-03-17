addpath(genpath('.'))   
allExpirements = dir(['expirements',filesep,'*.m']);

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
        clear config
        
        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

        [config,gpuFunc] = preprocessConfig(config);
        dimVec = [numel(config.focalPointsL.base), ...
            numel(config.focalPointsL.base), ...
            numel(config.focalPointsL.plain), ...
            numel(config.focalDirectionsL.base), ...
            numel(config.focalDirectionsL.base), ...
            numel(config.focalDirectionsV.base), ...
            numel(config.focalDirectionsV.base)
            ];
        tic
        [u,us] = run_rendering_mog(config,gpuFunc);
        u = reshape(u,dimVec);
        us = reshape(us,dimVec);
        t = toc

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'u','us','config','t');
    end
end
