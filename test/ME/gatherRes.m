clear

%% vmf
data.dirName = 'res_vmf';

expirementsRes = struct;

resFiles = dir([data.dirName,filesep,'*.mat']);


% randPixelWeight1e3 = 1e3 / (1e3 * (size(resFiles,1)/2));

for fileNum = 1:1:size(resFiles,1)
    load([data.dirName,filesep,resFiles(fileNum).name]);
    disp([num2str(fileNum),'/',num2str(size(resFiles,1))]);
    projectName = config.projectName;
    
    
    if(~isfield(expirementsRes,projectName))
         newFieldStruct.C = 0;
         newFieldStruct.Cs = 0;
         newFieldStruct.config = config;
         newFieldStruct.itersNum = 0;
         newFieldStruct.xRep = 0;
         
         newFieldStruct.fileCount = 0;
         
         expirementsRes = setfield(expirementsRes,projectName,newFieldStruct);
    end
    
    currField = getfield(expirementsRes,projectName);
    currField.C = currField.C + C;
    currField.itersNum = currField.itersNum + numOfIterations;
    currField.Cs = currField.Cs + Cs;
    currField.fileCount = currField.fileCount + 1;
    currField.xRep = currField.xRep + xRep;
    expirementsRes = setfield(expirementsRes,projectName,currField);
    
end


vmf = expirementsRes;
clear expirementsRes;
clear newFieldStruct;

%% gaussian

data.dirName = 'res_gaussian';

expirementsRes = struct;

resFiles = dir([data.dirName,filesep,'*.mat']);


% randPixelWeight1e3 = 1e3 / (1e3 * (size(resFiles,1)/2));

for fileNum = 1:1:size(resFiles,1)
    load([data.dirName,filesep,resFiles(fileNum).name]);
    disp([num2str(fileNum),'/',num2str(size(resFiles,1))]);
    projectName = config.projectName;
    
    
    if(~isfield(expirementsRes,projectName))
         newFieldStruct.C = 0;
         newFieldStruct.Cs = 0;
         newFieldStruct.config = config;
         newFieldStruct.itersNum = 0;
         
         newFieldStruct.fileCount = 0;
         
         expirementsRes = setfield(expirementsRes,projectName,newFieldStruct);
    end
    
    currField = getfield(expirementsRes,projectName);
    currField.C = currField.C + C;
    currField.itersNum = currField.itersNum + numOfIterations;
    currField.Cs = currField.Cs + Cs;
    currField.fileCount = currField.fileCount + 1;
    expirementsRes = setfield(expirementsRes,projectName,currField);
    
end

gaussian = expirementsRes;

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
         newFieldStruct.C = 0;
         newFieldStruct.Cs = 0;
         newFieldStruct.config = config;
         newFieldStruct.itersNum = 0;
         
         newFieldStruct.fileCount = 0;
         
         expirementsRes = setfield(expirementsRes,projectName,newFieldStruct);
    end
    
    currField = getfield(expirementsRes,projectName);
    currField.C = currField.C + C;
    currField.itersNum = currField.itersNum + numOfIterations;
    currField.Cs = currField.Cs + Cs;
    currField.fileCount = currField.fileCount + 1;
    expirementsRes = setfield(expirementsRes,projectName,currField);
    
end

ff = expirementsRes;

save('res.mat','gaussian','vmf','ff','-v7.3');