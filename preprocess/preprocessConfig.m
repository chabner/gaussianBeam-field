function [config] = preprocessConfig(config)
%% Sample
config.box_min = [-config.boxAxial/2;-config.boxAxial/2;-config.boxDepth/2];
config.box_max = [config.boxAxial/2;config.boxAxial/2;config.boxDepth/2];

if(isfield(config,'boxShift'))
    config.box_min = config.box_min + config.boxShift;
    config.box_max = config.box_max + config.boxShift;
end

if(~isfield(config,'multiplePaths'))
    config.multiplePaths = 1;
end

correlationActive = isfield(config.nf.focalPointsL,'x_2');

%% Parameters
currDim = 4;

paramsList = fieldnames(config.nf.parameters);

for paramNum = 1:1:numel(paramsList)
    param = config.nf.parameters.(paramsList{paramNum});
    param = param(:);
    permuteVec = 1:1:currDim;
    permuteVec(1) = currDim; permuteVec(end) = 1;
    param = permute(param,permuteVec);
    
    config.nf.parameters.(paramsList{paramNum}) = param;
    currDim = currDim + 1;
end

if(isfield(config,'parameters'))
    paramsList = fieldnames(config.parameters);

    for paramNum = 1:1:numel(paramsList)
        param = config.parameters.(paramsList{paramNum});
        param = param(:);
        permuteVec = 1:1:currDim;
        permuteVec(1) = currDim; permuteVec(end) = 1;
        param = permute(param,permuteVec);

        config.parameters.(paramsList{paramNum}) = param;
        currDim = currDim + 1;
    end
end

% eval wavelength
if(isa(config.wavelenght,'function_handle'))
    config.evalWavelenght = eval(['config.wavelenght',parseFunc(config.wavelenght,'config.parameters.')]);
else
    config.evalWavelenght = config.wavelenght;
end

%% Aperture
% evaluate

config.nf.focalPointsL.eval.x = eval(['config.nf.focalPointsL.x',parseFunc(config.nf.focalPointsL.x)]);
config.nf.focalPointsL.eval.y = eval(['config.nf.focalPointsL.y',parseFunc(config.nf.focalPointsL.y)]);
config.nf.focalPointsL.eval.z = eval(['config.nf.focalPointsL.z',parseFunc(config.nf.focalPointsL.z)]);

config.nf.focalPointsV.eval.x = eval(['config.nf.focalPointsV.x',parseFunc(config.nf.focalPointsV.x)]);
config.nf.focalPointsV.eval.y = eval(['config.nf.focalPointsV.y',parseFunc(config.nf.focalPointsV.y)]);
config.nf.focalPointsV.eval.z = eval(['config.nf.focalPointsV.z',parseFunc(config.nf.focalPointsV.z)]);

config.nf.focalDirectionsL.eval.x = eval(['config.nf.focalDirectionsL.x',parseFunc(config.nf.focalDirectionsL.x)]);
config.nf.focalDirectionsL.eval.y = eval(['config.nf.focalDirectionsL.y',parseFunc(config.nf.focalDirectionsL.y)]);
config.nf.focalDirectionsL.eval.z = eval(['config.nf.focalDirectionsL.z',parseFunc(config.nf.focalDirectionsL.z)]);

config.nf.focalDirectionsV.eval.x = eval(['config.nf.focalDirectionsV.x',parseFunc(config.nf.focalDirectionsV.x)]);
config.nf.focalDirectionsV.eval.y = eval(['config.nf.focalDirectionsV.y',parseFunc(config.nf.focalDirectionsV.y)]);
config.nf.focalDirectionsV.eval.z = eval(['config.nf.focalDirectionsV.z',parseFunc(config.nf.focalDirectionsV.z)]);

config.apertureVmf_l = movmfAperture(config.evalWavelenght, config.nf.mask_varL,-1,...
    config.nf.focalPointsL.eval.x    ,config.nf.focalPointsL.eval.y    ,config.nf.focalPointsL.eval.z,...
    config.nf.focalDirectionsL.eval.x,config.nf.focalDirectionsL.eval.y,config.nf.focalDirectionsL.eval.z);

config.apertureVmf_v = movmfAperture(config.evalWavelenght, config.nf.mask_varV,1,...
    config.nf.focalPointsV.eval.x    ,config.nf.focalPointsV.eval.y    ,config.nf.focalPointsV.eval.z,...
    config.nf.focalDirectionsV.eval.x,config.nf.focalDirectionsV.eval.y,config.nf.focalDirectionsV.eval.z);

if(correlationActive)
    config.nf.focalPointsL.eval.x_2 = eval(['config.nf.focalPointsL.x_2',parseFunc(config.nf.focalPointsL.x_2)]);
    config.nf.focalPointsL.eval.y_2 = eval(['config.nf.focalPointsL.y_2',parseFunc(config.nf.focalPointsL.y_2)]);
    config.nf.focalPointsL.eval.z_2 = eval(['config.nf.focalPointsL.z_2',parseFunc(config.nf.focalPointsL.z_2)]);

    config.nf.focalPointsV.eval.x_2 = eval(['config.nf.focalPointsV.x_2',parseFunc(config.nf.focalPointsV.x_2)]);
    config.nf.focalPointsV.eval.y_2 = eval(['config.nf.focalPointsV.y_2',parseFunc(config.nf.focalPointsV.y_2)]);
    config.nf.focalPointsV.eval.z_2 = eval(['config.nf.focalPointsV.z_2',parseFunc(config.nf.focalPointsV.z_2)]);

    config.nf.focalDirectionsL.eval.x_2 = eval(['config.nf.focalDirectionsL.x_2',parseFunc(config.nf.focalDirectionsL.x_2)]);
    config.nf.focalDirectionsL.eval.y_2 = eval(['config.nf.focalDirectionsL.y_2',parseFunc(config.nf.focalDirectionsL.y_2)]);
    config.nf.focalDirectionsL.eval.z_2 = eval(['config.nf.focalDirectionsL.z_2',parseFunc(config.nf.focalDirectionsL.z_2)]);

    config.nf.focalDirectionsV.eval.x_2 = eval(['config.nf.focalDirectionsV.x_2',parseFunc(config.nf.focalDirectionsV.x_2)]);
    config.nf.focalDirectionsV.eval.y_2 = eval(['config.nf.focalDirectionsV.y_2',parseFunc(config.nf.focalDirectionsV.y_2)]);
    config.nf.focalDirectionsV.eval.z_2 = eval(['config.nf.focalDirectionsV.z_2',parseFunc(config.nf.focalDirectionsV.z_2)]);

    apertureVmf_l_2 = movmfAperture(config.evalWavelenght, config.nf.mask_varL,-1,...
        config.nf.focalPointsL.eval.x_2    ,config.nf.focalPointsL.eval.y_2    ,config.nf.focalPointsL.eval.z_2,...
        config.nf.focalDirectionsL.eval.x_2,config.nf.focalDirectionsL.eval.y_2,config.nf.focalDirectionsL.eval.z_2);

    apertureVmf_v_2 = movmfAperture(config.evalWavelenght, config.nf.mask_varV,1,...
        config.nf.focalPointsV.eval.x_2    ,config.nf.focalPointsV.eval.y_2    ,config.nf.focalPointsV.eval.z_2,...
        config.nf.focalDirectionsV.eval.x_2,config.nf.focalDirectionsV.eval.y_2,config.nf.focalDirectionsV.eval.z_2);
    
    config.apertureVmf_l = movmfUnite(config.apertureVmf_l,apertureVmf_l_2,2);
    config.apertureVmf_v = movmfUnite(config.apertureVmf_v,apertureVmf_v_2,2);
end

%% Check legal directions
if( ...
        ~isreal(config.apertureVmf_l.dir1) ||  ...
        ~isreal(config.apertureVmf_l.dir2) || ....
        ~isreal(config.apertureVmf_l.dir3) || ....
        ~isreal(config.apertureVmf_v.dir1) || ....
        ~isreal(config.apertureVmf_v.dir2) || ....
        ~isreal(config.apertureVmf_v.dir3) )
    error('refocused directions cannot be ilegal');
end

dirlAbs = config.apertureVmf_l.dir1.^2 + config.apertureVmf_l.dir2.^2 + config.apertureVmf_l.dir3.^2;
dirlAbs = dirlAbs(:);
if( any(dirlAbs < 0.9999999 | dirlAbs > 1.0000001))
    error('refocused directions cannot be ilegal');
end

dirvAbs = config.apertureVmf_v.dir1.^2 + config.apertureVmf_v.dir2.^2 + config.apertureVmf_v.dir3.^2;
dirvAbs = dirvAbs(:);
if( any(dirvAbs < 0.9999999 | dirvAbs > 1.0000001))
    error('refocused directions cannot be ilegal');
end

%% Scattering function - entry
entry.sctType = config.sctType;
entry.dim = config.dimNum;
entry.vmf_k = config.vmf_k;
entry.vmf_iterations = config.vmf_iterations;

if(config.sctType == 3)
    entry.g = config.g;
    entry.forwardWeight = config.forwardWeight;
    entry.vmf_samples = config.vmf_samples;
    vmf_samples = entry.vmf_samples;
end

if(config.sctType == 4)
    entry.g = config.g;
    entry.forwardWeight = config.forwardWeight;
    entry.vmfKappaG = config.vmfKappaG;
    vmf_samples = 0;
end

if(config.sctType == 2)
    entry.pdf = config.pdf;
    vmf_samples = 0;
end

%% Scattering function - throughput
entryHash = DataHash(entry);

preprocessPath = which('vmfCache.mat');
if(isempty(preprocessPath))
    if(config.sctType == 3)
        config.ampfunc.g = config.g;
        config.ampfunc.forwardWeight = config.forwardWeight;
    end
    
    if(config.sctType == 4)
        config.ampfunc.g = config.g;
        config.ampfunc.forwardWeight = config.forwardWeight;
        config.ampfunc.vmfKappaG = config.vmfKappaG;
    end

    if(config.sctType == 2)
        theta = linspace(0,pi,numel(config.pdf));
        f = config.pdf;
        
        config.ampfunc.cs = mean((abs(f).^2).*sin(theta))*2*pi^2;
        config.ampfunc.samplePdf = (abs(f(:).').^2) .* sin(theta) ./ ...
            sum((abs(f(:).').^2) .* sin(theta));
        config.ampfunc.sampleCdf = cumsum(config.ampfunc.samplePdf);
        
        config.ampfunc.evalPdf = (abs(f(:).').^2) ./ config.ampfunc.cs;
        config.ampfunc.evalAmp = config.ampfunc.evalPdf .^ 0.5;
        
        config.ampfunc.sampleIcdf = invCDFvec(config.ampfunc.sampleCdf);
    end
    
    config.movmf = movmfBuild(config.sctType, config.ampfunc, config.dimNum, config.vmf_k, vmf_samples, config.vmf_iterations, config.useGpu);
    
    vmfCache(1).entryHash = entryHash;
    vmfCache(1).entry = entry;
    vmfCache(1).ampfunc = config.ampfunc;
    vmfCache(1).movmf = config.movmf;
    
    config.cacheIdx = 1;
    
    currentFilePath = which('preprocessConfig.m');
    currentFileFolder = regexp(currentFilePath,'C*[\/\\]preprocessConfig.m');
    
    save([currentFilePath(1:currentFileFolder),'vmfCache.mat'],'vmfCache');
else
    load(preprocessPath,'vmfCache')
    
    isInCache = false;
    for idxNum = 1:1:numel(vmfCache)
        if(strcmp(entryHash,vmfCache(idxNum).entryHash))
            config.movmf = vmfCache(idxNum).movmf;
            config.ampfunc = vmfCache(idxNum).ampfunc;
            config.cacheIdx = idxNum;
            isInCache = true;
            break;
        end
    end
    
    if(~isInCache)
        lastElem = numel(vmfCache);
        config.cacheIdx = lastElem + 1;
        
        if(config.sctType == 3)
            config.ampfunc.g = config.g;
            config.ampfunc.forwardWeight = config.forwardWeight;
        end
        
        if(config.sctType == 4)
            config.ampfunc.g = config.g;
            config.ampfunc.forwardWeight = config.forwardWeight;
            config.ampfunc.vmfKappaG = config.vmfKappaG;
        end

        if(config.sctType == 2)
            theta = linspace(0,pi,numel(config.pdf));
            f = config.pdf;

            config.ampfunc.cs = mean((abs(f).^2).*sin(theta))*2*pi^2;
            config.ampfunc.samplePdf = (abs(f(:).').^2) .* sin(theta) ./ ...
                sum((abs(f(:).').^2) .* sin(theta));
            config.ampfunc.sampleCdf = cumsum(config.ampfunc.samplePdf);

            config.ampfunc.evalPdf = (abs(f(:).').^2) ./ config.ampfunc.cs;
            config.ampfunc.evalAmp = config.ampfunc.evalPdf .^ 0.5;

            config.ampfunc.sampleIcdf = invCDFvec(config.ampfunc.sampleCdf);
        end
        
        config.movmf = movmfBuild(config.sctType, config.ampfunc, config.dimNum, config.vmf_k, vmf_samples, config.vmf_iterations, config.useGpu);
        
        vmfCache(lastElem+1).entryHash = entryHash;
        vmfCache(lastElem+1).entry = entry;
        vmfCache(lastElem+1).ampfunc = config.ampfunc;
        vmfCache(lastElem+1).movmf = config.movmf;
    
        save(preprocessPath,'vmfCache');
    end
end

%% Importance sampling

if(mod(config.sampleFlag,10) == 7)
    return;
end

% Sample 100 different samples, and take the median of px, in order to
% estimate bad samples
% pxItersNum = 1e1;
% pxItersNum = 1e3;
pxItersNum = 1e0;
minTorr = 1e-3;

box_w = config.box_max-config.box_min;
V=prod(box_w);

% if(mod(config.sampleFlag,10) == 3)
%     config.smpPreprocess = preprocess_smpVmfBeam(config,lightNum,1e4);
% end

if(mod(config.sampleFlag,10) == 4)
    config.smpPreprocess = preprocess_smpVmfBeamSum(config,config.apertureVmf_l);
end

switch mod(config.sampleFlag,10)
    case 1
        % uniform distribution
        px = 1;
    case 2
        % exponential distribution
        [~,px]=expSmpX(config.box_min,config.box_max,[0;0;1],1/config.MFP,pxItersNum);
    case 3
%         % gaussian sampling
%         apertureVmf_l = movmfAperture(config.mask_varL,config.focalPointsL{lightNum},-1,config.focalDirectionsL{lightNum});
%         [~,px(iterNum)] = smpVmfBeam(apertureVmf_l,config.smpPreprocess{lightNum},config.box_min,config.box_max);
%         px(iterNum) = px(iterNum) * V;
    case 4
        % gaussian sampling
        [~,px] = smpVmfBeamSum(config.apertureVmf_l,config.smpPreprocess,config.box_min,config.box_max,pxItersNum);
        px = px * V;
    case 6
        % known scattering positions
        px = 1;
end


config.smpPreprocess(1).pxMin = median(px,6) * minTorr;

end
