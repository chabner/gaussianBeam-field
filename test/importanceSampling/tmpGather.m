data.dirName = 'res';

expirementsRes1e3 = struct;
expirementsRes1e4 = struct;
expirementsRes1e5 = struct;
expirementsRes = struct;
expirementsResRandPixel1e5 = struct;
expirementsResRandPixel1e4 = struct;
expirementsResRandPixel1e3 = struct;
lastRound = struct;

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
         
         newFieldStruct.randPixelWeight1e5 = 1e5 / (1e3 * newFieldStruct.filesNum);
         newFieldStruct.randPixelWeight1e4 = 1e4 / (1e3 * newFieldStruct.filesNum);
         
         expirementsRes = setfield(expirementsRes,projectName,newFieldStruct);
         expirementsRes1e3 = setfield(expirementsRes1e3,projectName,newFieldStruct);
         expirementsRes1e4 = setfield(expirementsRes1e4,projectName,newFieldStruct);
         expirementsRes1e5 = setfield(expirementsRes1e5,projectName,newFieldStruct);
         expirementsResRandPixel1e5 = setfield(expirementsResRandPixel1e5,projectName,newFieldStruct);
         expirementsResRandPixel1e4 = setfield(expirementsResRandPixel1e4,projectName,newFieldStruct);
         expirementsResRandPixel1e3 = setfield(expirementsResRandPixel1e3,projectName,newFieldStruct);
         lastRound = setfield(lastRound,projectName,newFieldStruct);
    end
    
    currField = getfield(expirementsRes,projectName);
    currField.C = currField.C + C;
    currField.itersNum = currField.itersNum + numOfIterations;
    currField.Cs = currField.Cs + Cs;
    currField.fileCount = currField.fileCount + 1;
    currField.xRep = currField.xRep + xRep;
    expirementsRes = setfield(expirementsRes,projectName,currField);
    
    currField = getfield(expirementsRes1e3,projectName);
    if(currField.itersNum < 1e3)
        currField.C = currField.C + C;
        currField.itersNum = currField.itersNum + numOfIterations;
        currField.Cs = currField.Cs + Cs;
        currField.fileCount = currField.fileCount + 1;
        currField.xRep = currField.xRep + xRep;
        expirementsRes1e3 = setfield(expirementsRes1e3,projectName,currField);
    end
    
    currField = getfield(expirementsRes1e4,projectName);
    if(currField.itersNum < 1e4)
        currField.C = currField.C + C;
        currField.itersNum = currField.itersNum + numOfIterations;
        currField.Cs = currField.Cs + Cs;
        currField.fileCount = currField.fileCount + 1;
        currField.xRep = currField.xRep + xRep;
        expirementsRes1e4 = setfield(expirementsRes1e4,projectName,currField);
    end
    
    currField = getfield(expirementsRes1e5,projectName);
    if(currField.itersNum < 1e5)
        currField.C = currField.C + C;
        currField.itersNum = currField.itersNum + numOfIterations;
        currField.Cs = currField.Cs + Cs;
        currField.fileCount = currField.fileCount + 1;
        currField.xRep = currField.xRep + xRep;
        expirementsRes1e5 = setfield(expirementsRes1e5,projectName,currField);
    end
    
    currField = getfield(expirementsResRandPixel1e5,projectName);
    randPixels = double(rand(size(C)) < currField.randPixelWeight1e5);
    
    currField.C = currField.C + randPixels .* C;
    currField.itersNum = currField.itersNum + randPixels * numOfIterations;
    currField.Cs = currField.Cs + randPixels .* Cs;
    currField.fileCount = currField.fileCount + 1;
    currField.xRep = currField.xRep + randPixels .* xRep;
%     xRep = xRep(:);
%     currField.xRep = currField.xRep + randPixels .* permute(xRep,[3,2,1]);
    expirementsResRandPixel1e5 = setfield(expirementsResRandPixel1e5,projectName,currField);
    
    currField = getfield(expirementsResRandPixel1e4,projectName);
    randPixels = double(rand(size(C)) < currField.randPixelWeight1e4);
    
    currField.C = currField.C + randPixels .* C;
    currField.itersNum = currField.itersNum + randPixels * numOfIterations;
    currField.Cs = currField.Cs + randPixels .* Cs;
    currField.fileCount = currField.fileCount + 1;
    currField.xRep = currField.xRep + randPixels .* xRep;
%     xRep = xRep(:);
%     currField.xRep = currField.xRep + randPixels .* permute(xRep,[3,2,1]);
    expirementsResRandPixel1e4 = setfield(expirementsResRandPixel1e4,projectName,currField);
    
    currField = getfield(expirementsResRandPixel1e3,projectName);
    currField.fileCount = currField.fileCount + 1;
    
    randPixels = (currField.pixelFileNum == currField.fileCount);
    
    currField.C = currField.C + randPixels .* C;
    currField.itersNum = currField.itersNum + randPixels * numOfIterations;
    currField.Cs = currField.Cs + randPixels .* Cs;
%     xRep = xRep(:);
%     currField.xRep = currField.xRep + randPixels .* permute(xRep,[3,2,1]);
    expirementsResRandPixel1e3 = setfield(expirementsResRandPixel1e3,projectName,currField);
    
    currField = getfield(lastRound,projectName);    
    currField.C = C;
    currField.itersNum = 1e3;
    currField.Cs = Cs;
%     xRep = xRep(:);
%     currField.xRep = currField.xRep + randPixels .* permute(xRep,[3,2,1]);
    lastRound = setfield(lastRound,projectName,currField);
    
end

% % fill from the last round for empty cells
% fn = fieldnames(expirementsResRandPixel1e3);
% for k=1:numel(fn)
%     res = expirementsResRandPixel1e3.(fn{k});
%     lr = lastRound.(fn{k});
%     emptyDots = (res.itersNum == 0);
%     
%     res.C(emptyDots) = lr.C(emptyDots);
%     res.itersNum(emptyDots) = lr.itersNum;
%     res.I_1(emptyDots) = lr.I_1(emptyDots);
%     res.I_2(emptyDots) = lr.I_2(emptyDots);
%     res.Cs(emptyDots) = lr.Cs(emptyDots);
%     res.Is_1(emptyDots) = lr.Is_1(emptyDots);
%     res.Is_2(emptyDots) = lr.Is_2(emptyDots);
% %     xRep = xRep(:);
% %     currField.xRep = currField.xRep + randPixels .* permute(xRep,[3,2,1]);
%     expirementsResRandPixel1e3.(fn{k}) = res;
% end


save([data.dirName,'_tmp.mat'],'expirementsRes','expirementsRes1e3','expirementsRes1e4', ...
    'expirementsRes1e5','expirementsResRandPixel1e3','expirementsResRandPixel1e4',...
    'expirementsResRandPixel1e5','-v7.3');