clear

allExpirements = dir('config*.m');

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
            clear config
            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

            config = preprocessConfig(config);
            
            C = 0;
            Cs = 0;
            xRep = 0;
            
            tic
            for corrIter = 1:1:1e3
                [u,us,current_xRep] = run_rendering(config);
                C = C + u(:,:,1) .* conj(u);
                Cs = Cs + us(:,:,1) .* conj(us);
                xRep = xRep + current_xRep;
            end
            t = toc
            numOfIterations = 1e3;

            T = datetime('now','Format','ddMM_HHmmss_SSS');
            save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
                'C','Cs','config','t','xRep','numOfIterations');
    end
end
