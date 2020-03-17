clear

%% nf
data.dirName = 'res_nf';

expirementsRes = struct;

resFiles = dir([data.dirName,filesep,'*.mat']);


% randPixelWeight1e3 = 1e3 / (1e3 * (size(resFiles,1)/2));

for fileNum = 1:1:size(resFiles,1)
    load([data.dirName,filesep,resFiles(fileNum).name]);
    disp([num2str(fileNum),'/',num2str(size(resFiles,1))]);
    projectName = config.projectName;
    
    
    if(~isfield(expirementsRes,projectName))
         newFieldStruct.C_nf = 0;
         newFieldStruct.Cs_nf = 0;
         newFieldStruct.config = config;
         newFieldStruct.itersNum = 0;
         newFieldStruct.xRep = 0;
         
         newFieldStruct.fileCount = 0;
         
         expirementsRes = setfield(expirementsRes,projectName,newFieldStruct);
    end
    
    currField = getfield(expirementsRes,projectName);
    currField.C_nf = currField.C_nf + C_nf;
    currField.itersNum = currField.itersNum + numOfIterations;
    currField.Cs_nf = currField.Cs_nf + Cs_nf;
    currField.fileCount = currField.fileCount + 1;
    currField.xRep = currField.xRep + xRep;
    expirementsRes = setfield(expirementsRes,projectName,currField);
    
end


nf = expirementsRes;
clear expirementsRes;
clear newFieldStruct;

%% ff

data.dirName = 'res_ff';

expirementsRes = struct;

resFiles = dir([data.dirName,filesep,'*.mat']);


% randPixelWeight1e3 = 1e3 / (1e3 * (size(resFiles,1)/2));

for fileNum = 1:1:size(resFiles,1)
    load([data.dirName,filesep,resFiles(fileNum).name]);
    disp([num2str(fileNum),'/',num2str(size(resFiles,1))]);
    projectName = config.projectName;
    
    
    if(~isfield(expirementsRes,projectName))
         newFieldStruct.C_ff = 0;
         newFieldStruct.Cs_ff = 0;
         newFieldStruct.config = config;
         newFieldStruct.itersNum = 0;
         
         newFieldStruct.fileCount = 0;
         
         expirementsRes = setfield(expirementsRes,projectName,newFieldStruct);
    end
    
    currField = getfield(expirementsRes,projectName);
    currField.C_ff = currField.C_ff + C_ff;
    currField.itersNum = currField.itersNum + numOfIterations;
    currField.Cs_ff = currField.Cs_ff + Cs_ff;
    currField.fileCount = currField.fileCount + 1;
    expirementsRes = setfield(expirementsRes,projectName,currField);
    
end

ff = expirementsRes;

save('res.mat','ff','nf','l_base','-v7.3');