clear

allExpirements = dir('config*.m');

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
%         %% ff
%             clear config
% 
%             run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);
%             
%             D = abs(config.focalPointsL.plain) + abs(config.boxDepth/2);
%             v_base = atan(config.focalPointsV.base / D);
%             l_base = atan(lBase / D);
%             config.sampleFlag = 32;
%             
%             [config] = preprocessConfig_ff(config);
%             
            numOfIterations = 1e4;
%             
%             tic
%             [C_ff,Cs_ff] = farField_correlation(config,l_base,v_base,numOfIterations);
%             t = toc
% 
%             T = datetime('now','Format','ddMM_HHmmss_SSS');
%             save(['res_ff',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
%                 'C_ff','Cs_ff','config','t','numOfIterations','l_base','v_base');
      %% nf
            clear config

            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);
            config.sampleFlag = 42;
            
            [config,gpuFunc] = preprocessConfig(config);
            
            C_nf = 0;
            Cs_nf = 0;
            xRep = 0;
            
            Nl = numel(config.focalPointsL.base);
            Nv = numel(config.focalPointsV.base);
            
            tic
            for corrIter = 1:1:numOfIterations
                [u,us,current_xRep] = run_rendering(config,gpuFunc);
                u = reshape(u,[Nl,Nv,Nv]); u = permute(u,[2,3,1]);
                us = reshape(us,[Nl,Nv,Nv]); us = permute(us,[2,3,1]);
                
                C_nf = C_nf + u(:,:,(1:Nl/2)) .* conj(u(:,:,((1+Nl/2):end)));
                Cs_nf = Cs_nf + us(:,:,(1:Nl/2)) .* conj(us(:,:,((1+Nl/2):end)));
                xRep = xRep + current_xRep;
            end
            t = toc
            
            C_nf = gather(C_nf);
            Cs_nf = gather(Cs_nf);

            T = datetime('now','Format','ddMM_HHmmss_SSS');
            save(['res_nf',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
                'C_nf','Cs_nf','config','t','xRep','numOfIterations');
    end
end
