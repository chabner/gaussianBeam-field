   
allExpirements = dir(['expirements',filesep,'*.m']);

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
        clear config
        
        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

        config = preprocessConfig_onGrid(config);
        dimVec = [numel(config.focalPoints_base), ...
            numel(config.focalPoints_base), ...
            numel(config.focalPoints_plain), ...
            numel(config.focalLDirections), ...
            numel(config.focalLDirections), ...
            numel(config.focalVDirections), ...
            numel(config.focalVDirections)
            ];
        tic
        [u,us] = run_rendering(config);
        u = reshape(u,dimVec);
        us = reshape(us,dimVec);
        t = toc

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'u','us','config','t');
    end
end
