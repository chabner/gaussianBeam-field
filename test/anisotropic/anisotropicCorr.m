clear

allExpirements = dir('config*.m');

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
        %% vMF
        
        clear config

        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

        [config,gpuFunc] = preprocessConfig(config);

        C = 0;
        Cs = 0;
        xRep = 0;

        Nl = numel(config.focalPointsL.base);
        Nv = numel(config.focalPointsV.base);
        
        numOfIterations = 1e4;
        tic
        for corrIter = 1:1:numOfIterations
            [u,us,current_xRep] = run_rendering(config,gpuFunc);
            u = reshape(u,[Nl,Nv,Nv]); u = permute(u,[2,3,1]);
            us = reshape(us,[Nl,Nv,Nv]); us = permute(us,[2,3,1]);

            C = C + u(:,:,(1:Nl/2)) .* conj(u(:,:,((1+Nl/2):end)));
            Cs = Cs + us(:,:,(1:Nl/2)) .* conj(us(:,:,((1+Nl/2):end)));
            xRep = xRep + current_xRep;
        end
        t = toc

        C = gather(C);
        Cs = gather(Cs);

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res_vmf',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'C','Cs','config','t','xRep','numOfIterations');
        
        clear config
        
        %% Gaussian

        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

        [config,gpuFunc] = preprocessConfig(config);

        C = 0;
        Cs = 0;
        xRep = 0;

        Nl = numel(config.focalPointsL.base);
        Nv = numel(config.focalPointsV.base);
        
        numOfIterations = 1e4;

        tic
        for corrIter = 1:1:numOfIterations
            [u,us,current_xRep] = run_rendering_mog(config,gpuFunc);
            u = reshape(u,[Nl,Nv,Nv]); u = permute(u,[2,3,1]);
            us = reshape(us,[Nl,Nv,Nv]); us = permute(us,[2,3,1]);

            C = C + u(:,:,(1:Nl/2)) .* conj(u(:,:,((1+Nl/2):end)));
            Cs = Cs + us(:,:,(1:Nl/2)) .* conj(us(:,:,((1+Nl/2):end)));
            xRep = xRep + current_xRep;
        end
        t = toc

        C = gather(C);
        Cs = gather(Cs);

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res_gaussian',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'C','Cs','config','t','xRep','numOfIterations');
        
        clear config

        %% ff

        run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);
        
        config.g0 = config.g;
        config.sampleFlag = 32;

        config = preprocessConfig_ff(config);
        
        Nl = numel(config.focalPointsL.base);
        Nv = numel(config.focalPointsV.base);
        
        l_base = 2*abs(config.focalPointsL.base(:,1:Nl/2)) ./ config.boxDepth;
        v_base = config.focalPointsV.base ./ config.boxDepth;
        numOfIterations = 1e4;
        
        tic
        [C,Cs] = farField_correlation(config,l_base,v_base,numOfIterations) ;
        t = toc


        C = gather(C);
        Cs = gather(Cs);

        T = datetime('now','Format','ddMM_HHmmss_SSS');
        save(['res_ff',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
            'C','Cs','config','t','numOfIterations');
    end
end
