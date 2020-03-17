clear;

configFileList = dir(['res_good',filesep,'config*']);
Cs = 0;
C = 0;
itersNum = 0;

for fileInConfigList = 1:1:numel(configFileList)
    load([configFileList(fileInConfigList).folder,filesep,configFileList(fileInConfigList).name]);
    Cs = Cs + Cs_ff;
    C = C + C_ff;
    itersNum = itersNum + numOfIterations;
end

Cs_ff = Cs;
C_ff = C;
numOfIterations = itersNum;

save('gatheredRes','Cs_ff','C_ff','config','numOfIterations','l_base','v_base')