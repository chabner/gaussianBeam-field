clear

allExpirements = dir('config*.m');

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
            clear config
            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

            [config,gpuFunc] = preprocessConfig(config);
            
            C = 0;
            Cs = 0;
            xRep = 0;
            
            Nl = numel(config.focalPointsL.base);
            Nv = numel(config.focalPointsV.base);
            
            tic
            for corrIter = 1:1:1e4
                [u,us,current_xRep] = run_rendering(config,gpuFunc);
                u = reshape(u,[Nl,Nv,Nv]); u = permute(u,[2,3,1]);
                us = reshape(us,[Nl,Nv,Nv]); us = permute(us,[2,3,1]);
                
                C = C + u(:,:,(1:Nl/2)) .* conj(u(:,:,((1+Nl/2):end)));
                Cs = Cs + us(:,:,(1:Nl/2)) .* conj(us(:,:,((1+Nl/2):end)));
                xRep = xRep + current_xRep;
            end
            t = toc
            numOfIterations = 1e4;
            
            C = gather(C);
            Cs = gather(Cs);

            T = datetime('now','Format','ddMM_HHmmss_SSS');
            save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
                'C','Cs','config','t','xRep','numOfIterations');
    end
end
