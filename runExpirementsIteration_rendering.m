   
allExpirements = dir(['expirements',filesep,'*.m']);

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

        config = preprocessConfig(config);
        tic
        [u,us] = run_rendering(config);
        t = toc

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'u','us','config','t');
    end
end
