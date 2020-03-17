clear

allExpirements = dir(['exp',filesep,'config*.m']);
numOfIterations = 1e4;

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
        %% vmf
            clear config
            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

            [config,gpuFunc] = preprocessConfig(config);
            
            C = 0;
            Cs = 0;
            xRep = 0;
            
            Nl = numel(config.focalPointsL.base);
            Nv = numel(config.focalPointsV.base);
            
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
            
        %% gaussian
            clear config
            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

            [config,gpuFunc] = preprocessConfig(config);
            
            C = 0;
            Cs = 0;
            xRep = 0;
            
            Nl = numel(config.focalPointsL.base);
            Nv = numel(config.focalPointsV.base);
            
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
            
%         %% ff
%             clear config
%             
%             if(~isempty(regexp(allExpirements(expirementFile).name,'_ff','once')))
%                 run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);
% 
%                 Nl = numel(config.focalPointsL.base);
% 
%                 D = abs(config.focalPointsL.plain) + abs(config.boxDepth/2);
%                 v_base = atan(config.focalPointsV.base / D);
%                 l_base = atan(lBase / D);
%                 config.sampleFlag = 32;
% 
%                 [config] = preprocessConfig_ff(config);
% 
%                 tic
%                 [C,Cs] = farField_correlation(config,l_base,v_base,numOfIterations);
%                 t = toc
%                 
%                 C = gather(C);
%                 Cs = gather(Cs);
% 
%                 T = datetime('now','Format','ddMM_HHmmss_SSS');
%                 save(['res_ff',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
%                     'C','Cs','config','t','numOfIterations','l_base','v_base'); 
%             end
    end
end
