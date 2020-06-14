function [config] = preprocessConfig2D(config)
%% Sample
config.box_min = [-config.boxAxial/2;-config.boxDepth/2];
config.box_max = [config.boxAxial/2;config.boxDepth/2];

if(isfield(config,'boxShift'))
    config.box_min = config.box_min + config.boxShift;
    config.box_max = config.box_max + config.boxShift;
end

%% Aperture
% Focal points - L

[focalPointsL_X,focalPointsL_Z] = ndgrid( ...
    config.focalPointsL.base, ...
    config.focalPointsL.plain);
config.focalPointsL.vector = [focalPointsL_X(:).'; focalPointsL_Z(:).'];

% Focal points - V
[focalPointsV_X,focalPointsV_Z] = ndgrid( ...
    config.focalPointsV.base, ...
    config.focalPointsV.plain);
config.focalPointsV.vector = [focalPointsV_X(:).'; focalPointsV_Z(:).'];

% V shift
if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift'))
    [vShift_X,vShift_Z] = ndgrid( ...
        config.focalPointsVshift.vShift, ...
        0);
    
    config.focalPointsVshift.vShiftVector = [vShift_X(:).'; vShift_Z(:).'];
end

% Focal directions - L
config.focalDirectionsL.vector = [ ...
    (sin(config.focalDirectionsL.theta(:)).') ; ...
    (cos(config.focalDirectionsL.theta(:)).') ];
    
% Focal directions - V
config.focalDirectionsV.vector = [ ...
    (sin(config.focalDirectionsV.theta(:)).') ; ...
    (cos(config.focalDirectionsV.theta(:)).') ];

config.apertureVmf_l = movmfAperture2D( ...
    config.mask_varL , ...
    config.focalPointsL.vector , ...
    -1 , ...
    config.focalDirectionsL.vector , ...
    config.focalPointsL.dim , ...
    config.focalDirectionsL.dim);

if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift'))
    config.apertureVmf_v = movmfAperture2D( ...
        config.mask_varV , ...
        config.focalPointsV.vector , ...
        1 , ...
        config.focalDirectionsV.vector , ...
        config.focalPointsV.dim , ...
        config.focalDirectionsV.dim, ...
        config.focalPointsVshift.vShiftVector, ...
        config.focalPointsVshift.dim);
else
    config.apertureVmf_v = movmfAperture2D( ...
        config.mask_varV , ...
        config.focalPointsV.vector , ...
        1 , ...
        config.focalDirectionsV.vector , ...
        config.focalPointsV.dim , ...
        config.focalDirectionsV.dim);
end

%% Aperture - correlation
if(isfield(config.focalPointsL,'base_2'))
    % Focal points - L
    [focalPointsL_X,focalPointsL_Z] = ndgrid( ...
        config.focalPointsL.base_2, ...
        config.focalPointsL.plain_2);
    config.focalPointsL.vector_2 = [focalPointsL_X(:).'; focalPointsL_Z(:).'];

    % Focal points - V
    [focalPointsV_X,focalPointsV_Z] = ndgrid( ...
        config.focalPointsV.base_2, ...
        config.focalPointsV.plain_2);
    config.focalPointsV.vector_2 = [focalPointsV_X(:).'; focalPointsV_Z(:).'];

    % V shift
    if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift_2'))
        [vShift_X,vShift_Z] = ndgrid( ...
            config.focalPointsVshift.vShift_2, ...
            0);

        config.focalPointsVshift.vShiftVector_2 = [vShift_X(:).'; vShift_Z(:).'];
    end

    % Focal directions - L
    config.focalDirectionsL.vector_2 = [ ...
        (sin(config.focalDirectionsL.theta_2(:)).') ; ...
        (cos(config.focalDirectionsL.theta_2(:)).') ];

    % Focal directions - V
    config.focalDirectionsV.vector_2 = [ ...
        (sin(config.focalDirectionsV.theta_2(:)).') ; ...
        (cos(config.focalDirectionsV.theta_2(:)).') ];

    config.apertureVmf_l_2 = movmfAperture2D( ...
        config.mask_varL , ...
        config.focalPointsL.vector_2 , ...
        -1 , ...
        config.focalDirectionsL.vector_2 , ...
        config.focalPointsL.dim , ...
        config.focalDirectionsL.dim);

    if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift_2'))
        config.apertureVmf_v_2 = movmfAperture2D( ...
            config.mask_varV , ...
            config.focalPointsV.vector_2 , ...
            1 , ...
            config.focalDirectionsV.vector_2 , ...
            config.focalPointsV.dim , ...
            config.focalDirectionsV.dim, ...
            config.focalPointsVshift.vShiftVector_2, ...
            config.focalPointsVshift.dim);
    else
        config.apertureVmf_v_2 = movmfAperture2D( ...
            config.mask_varV , ...
            config.focalPointsV.vector_2 , ...
            1 , ...
            config.focalDirectionsV.vector_2 , ...
            config.focalPointsV.dim , ...
            config.focalDirectionsV.dim);
    end
end

%% Scattering function - g0
config.ampfunc0 = inf;

%% Scattering function - entry

if(~isfield(config,'albedo'))
    config.albedo = 1;
end

entry.sctType = config.sctType;
entry.dim = config.dimNum;
entry.vmf_k = config.vmf_k;
entry.vmf_iterations = config.vmf_iterations;
entry.albedo = config.albedo;

if(config.sctType == 3)
    entry.g = config.g;
    entry.forwardWeight = config.forwardWeight;
    entry.vmf_samples = config.vmf_samples;
    vmf_samples = entry.vmf_samples;
end

if(config.sctType == 2)
    entry.pdf = config.pdf;
    vmf_samples = 0;
end
%% Importance sampling

if(isfield(config.focalPointsL,'base_2'))
    apertureVmf_l = movmfUnite2D(config.apertureVmf_l, config.apertureVmf_l_2);
else
    apertureVmf_l = config.apertureVmf_l;
end

% Sample 100 different samples, and take the median of px, in order to
% estimate bad samples
pxItersNum = 100;
minTorr = 1e-3;

box_w = config.box_max-config.box_min;
V=prod(box_w);

if(mod(config.sampleFlag,10) == 4)
    config.smpPreprocess = preprocess_smpVmfBeamSum2D(config,1e4,apertureVmf_l);
end

px = zeros(1,pxItersNum,class(config.apertureVmf_l.alpha));

for iterNum = 1:1:pxItersNum
    switch mod(config.sampleFlag,10)
    case 1
        % uniform distribution
        px(iterNum) = 1;
    case 2
        % exponential distribution
        [~,px(iterNum)]=expSmpX(config.box_min,config.box_max,[0;1],1/config.attMFP);
    case 4
        % gaussian sampling
        [~,px(iterNum)] = smpVmfBeamSum2D(apertureVmf_l,config.smpPreprocess,config.box_min,config.box_max);
        px(iterNum) = px(iterNum) * V;
    end
end

config.smpPreprocess(1).pxMin = median(px) * minTorr;

%% Scattering function - throughput

entryHash = DataHash(entry);

preprocessPath = which('vmfCache.mat');
if(isempty(preprocessPath))
    if(config.sctType == 3)
        config.ampfunc.g = config.g;
        config.ampfunc.forwardWeight = config.forwardWeight;
    end

    if(config.sctType == 2)
        f = config.pdf;
        
        config.ampfunc.cs = (mean(abs(f).^2)*2*pi) * config.albedo;
        config.ampfunc.samplePdf = (abs(f(:).').^2) ./ sum((abs(f(:).').^2));
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
            isInCache = true;
            break;
        end
    end
    
    if(~isInCache)
        lastElem = numel(vmfCache);
        
        if(config.sctType == 3)
            config.ampfunc.g = config.g;
            config.ampfunc.forwardWeight = config.forwardWeight;
        end

        if(config.sctType == 2)
            f = config.pdf;

            config.ampfunc.cs = mean(abs(f).^2)*2*pi * config.albedo;
            config.ampfunc.samplePdf = (abs(f(:).').^2) ./ sum((abs(f(:).').^2));
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

end