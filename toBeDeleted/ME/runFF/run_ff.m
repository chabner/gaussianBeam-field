clear

allExpirements = dir('config*.m');

rng('shuffle')
for iterNum = 1:1:1e3
    for expirementFile = 1:1:numel(allExpirements)
            clear config
            load('MEcor_Jan192020_od1.mat')

            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);
            l_base = [0:40:200,300:100:1000,2000:1000:4000]/(20000+100);
            
            v_max=0.16;
            v_stp=v_max/20;
            v_base = (-v_max):v_stp:v_max;

            [config] = preprocessConfig_ff(config);
            
            C = 0;
            Cs = 0;
            xRep = 0;
            
            Nl = numel(config.focalPointsL.base);
            Nv = numel(config.focalPointsV.base);
            
            numOfIterations = 1e4;
            
            tic
            [C_ff,Cs_ff] = farField_correlation(config,l_base,v_base,numOfIterations);
            t = toc

            T = datetime('now','Format','ddMM_HHmmss_SSS');
            save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
                'C_ff','Cs_ff','config','t','numOfIterations','l_base','v_base');
    end
end
