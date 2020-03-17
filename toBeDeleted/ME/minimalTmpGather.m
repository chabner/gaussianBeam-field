data.dirName = 'res';

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
         
         newFieldStruct.filesNum = sum(numel(dir([data.dirName,filesep,'*',projectName,'*.mat'])));
         newFieldStruct.pixelFileNum = randi(newFieldStruct.filesNum,size(C));
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


save([data.dirName,'_tmp.mat'],'expirementsRes','-v7.3');