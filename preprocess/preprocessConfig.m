function [config,gpuFunc] = preprocessConfig(config)
%% Sample
config.box_min = [-config.boxAxial/2;-config.boxAxial/2;-config.boxDepth/2];
config.box_max = [config.boxAxial/2;config.boxAxial/2;config.boxDepth/2];

if(isfield(config,'boxShift'))
    config.box_min = config.box_min + config.boxShift;
    config.box_max = config.box_max + config.boxShift;
end

if(config.mcGpuL)
    config.box_min = gpuArray(config.box_min);
    config.box_max = gpuArray(config.box_max);
end

%% Aperture
% Focal points - L

if(config.focalPointsL.xyGrid)
    lyVec = config.focalPointsL.base;
else
    lyVec = 0;
end

[focalPointsL_X,focalPointsL_Y,focalPointsL_Z] = ndgrid( ...
    config.focalPointsL.base, ...
    lyVec, ...
    config.focalPointsL.plain);
config.focalPointsL.vector = [focalPointsL_X(:).'; focalPointsL_Y(:).'; focalPointsL_Z(:).'];

% Focal points - V
if(config.focalPointsV.xyGrid)
    vyVec = config.focalPointsV.base;
else
    vyVec = 0;
end

[focalPointsV_X,focalPointsV_Y,focalPointsV_Z] = ndgrid( ...
    config.focalPointsV.base, ...
    vyVec, ...
    config.focalPointsV.plain);
config.focalPointsV.vector = [focalPointsV_X(:).'; focalPointsV_Y(:).'; focalPointsV_Z(:).'];

% Focal directions - L
if(config.focalDirectionsL.xyGrid)
    ldiryVec = sin(config.focalDirectionsL.base);
else
    ldiryVec = 0;
end

[focalDirectionsL_X,focalDirectionsL_Y] = ndgrid( ...
    sin(config.focalDirectionsL.base), ...
    ldiryVec);

config.focalDirectionsL.vector = [ ...
    focalDirectionsL_X(:).' ; ...
    focalDirectionsL_Y(:).' ; ...
    sqrt(1 - (focalDirectionsL_X(:).').^2 - (focalDirectionsL_Y(:).').^2)];

% Focal directions - V
if(config.focalDirectionsV.xyGrid)
    vdiryVec = sin(config.focalDirectionsV.base);
else
    vdiryVec = 0;
end

[focalDirectionsV_X,focalDirectionsV_Y] = ndgrid( ...
    sin(config.focalDirectionsV.base), ...
    vdiryVec);

config.focalDirectionsV.vector = [ ...
    focalDirectionsV_X(:).' ; ...
    focalDirectionsV_Y(:).' ; ...
    sqrt(1 - (focalDirectionsV_X(:).').^2 - (focalDirectionsV_Y(:).').^2)];

% to gpu

if(config.mcGpuL)
    config.focalPointsL.vector = gpuArray(config.focalPointsL.vector);
    config.focalDirectionsL.vector = gpuArray(config.focalDirectionsL.vector);
end

if(config.mcGpuV)
    config.focalPointsV.vector = gpuArray(config.focalPointsV.vector);
    config.focalDirectionsV.vector = gpuArray(config.focalDirectionsV.vector);
end

config.apertureVmf_l = movmfAperture( ...
    config.mask_varL , ...
    config.focalPointsL.vector , ...
    -1 , ...
    config.focalDirectionsL.vector , ...
    config.focalPointsL.dim , ...
    config.focalDirectionsL.dim);

config.apertureVmf_v = movmfAperture( ...
    config.mask_varV , ...
    config.focalPointsV.vector , ...
    1 , ...
    config.focalDirectionsV.vector , ...
    config.focalPointsV.dim , ...
    config.focalDirectionsV.dim);

%% Scattering function
if(config.sctType == 3)
    config.sct_type = 3;
    config.ampfunc.g = config.g;
    config.ampfunc.forwardWeight = config.forwardWeight;
    config.ampfunc0 = inf;
end

%% Scattering function - throughput
preprocessPath = which('vmfCache.mat');
if(isempty(preprocessPath))
    config.movmf = movmfBuild(config.sctType, config.ampfunc, config.dimNum, config.vmf_k, config.vmf_samples, config.vmf_iterations, config.useGpu);
    
    vmfCache(1).g = config.g;
    vmfCache(1).forwardWeight = config.forwardWeight;
    vmfCache(1).vmf_k = config.vmf_k;
    vmfCache(1).vmf_iterations = config.vmf_iterations;
    vmfCache(1).vmf_samples = config.vmf_samples;
    vmfCache(1).movmf = config.movmf;
    
    currentFilePath = which('preprocessConfig.m');
    currentFileFolder = regexp(currentFilePath,'C*[\/\\]preprocessConfig.m');
    
    save([currentFilePath(1:currentFileFolder),'vmfCache.mat'],'vmfCache');
else
    load(preprocessPath,'vmfCache')
    
    isInCache = false;
    for idxNum = 1:1:numel(vmfCache)
        if( ...
            vmfCache(idxNum).g == config.g                             && ...
            vmfCache(idxNum).forwardWeight == config.forwardWeight     && ...
            vmfCache(idxNum).vmf_k == config.vmf_k                     && ...
            vmfCache(idxNum).vmf_iterations == config.vmf_iterations   && ...
            vmfCache(idxNum).vmf_samples == config.vmf_samples            ...
        )
            config.movmf = vmfCache(idxNum).movmf;
            isInCache = true;
            break;
        end
    end
    
    if(~isInCache)
        config.movmf = movmfBuild(config.sctType, config.ampfunc, config.dimNum, config.vmf_k, config.vmf_samples, config.vmf_iterations, config.useGpu);
        lastElem = numel(vmfCache);
    
        vmfCache(lastElem+1).g = config.g;
        vmfCache(lastElem+1).forwardWeight = config.forwardWeight;
        vmfCache(lastElem+1).vmf_k = config.vmf_k;
        vmfCache(lastElem+1).vmf_iterations = config.vmf_iterations;
        vmfCache(lastElem+1).vmf_samples = config.vmf_samples;
        vmfCache(lastElem+1).movmf = config.movmf;
    
        save(preprocessPath,'vmfCache');
    end
end

%% Importance sampling
% Sample 100 different samples, and take the median of px, in order to
% estimate bad samples
pxItersNum = 100;
minTorr = 1e-3;

box_w = config.box_max-config.box_min;
V=prod(box_w);

if(mod(config.sampleFlag,10) == 3)
    config.smpPreprocess = preprocess_smpVmfBeam(config,lightNum,1e4);
end

if(mod(config.sampleFlag,10) == 4)
    config.smpPreprocess = preprocess_smpVmfBeamSum(config,1e4);
end

px = zeros(1,pxItersNum,class(config.apertureVmf_l.alpha));

for iterNum = 1:1:pxItersNum
    switch mod(config.sampleFlag,10)
    case 1
        % uniform distribution
        px(iterNum) = 1;
    case 2
        % exponential distribution
        [~,px(iterNum)]=expSmpX(config.box_min,config.box_max,[0;0;1],config.attMFP);
    case 3
        % gaussian sampling
        apertureVmf_l = movmfAperture(config.mask_varL,config.focalPointsL{lightNum},-1,config.focalDirectionsL{lightNum});
        [~,px(iterNum)] = smpVmfBeam(apertureVmf_l,config.smpPreprocess{lightNum},config.box_min,config.box_max);
        px(iterNum) = px(iterNum) * V;
    case 4
        % gaussian sampling
        [~,px(iterNum)] = smpVmfBeamSum(config.apertureVmf_l,config.smpPreprocess,config.box_min,config.box_max);
        px(iterNum) = px(iterNum) * V;
    end
end

config.smpPreprocess(1).pxMin = median(px) * minTorr;

%% GPU code
if(config.mcGpuL || config.mcGpuV)
    gpuFunc.active = true;
    u_size = max(config.apertureVmf_l.dim,config.apertureVmf_v.dim);
    u_size = u_size(2:end);
    
    gpuFunc.integrateMult = parallel.gpu.CUDAKernel('integrateMult.ptx','integrateMult.cu');
    gpuFunc.integrateMult.GridSize = [ceil(prod(u_size)/gpuFunc.integrateMult.MaxThreadsPerBlock) 1 1];
    gpuFunc.integrateMult.ThreadBlockSize = [gpuFunc.integrateMult.MaxThreadsPerBlock 1 1];
    setConstantMemory(gpuFunc.integrateMult,'uDimProd',int32(cumprod(u_size)));
    
    lDim = config.apertureVmf_l.dim;
    lDim(1) = config.movmf.dim(1);
    setConstantMemory(gpuFunc.integrateMult,'lDim',int32(lDim));
    setConstantMemory(gpuFunc.integrateMult,'vDim',int32(config.apertureVmf_v.dim));
    setConstantMemory(gpuFunc.integrateMult,'lMixtureAlpha',config.movmf.alpha);
    
    gpuFunc.MSintegrateMult = parallel.gpu.CUDAKernel('MSintegrateMult.ptx','MSintegrateMult.cu');
    gpuFunc.MSintegrateMult.GridSize = [ceil(prod(u_size)/gpuFunc.MSintegrateMult.MaxThreadsPerBlock) 1 1];
    gpuFunc.MSintegrateMult.ThreadBlockSize = [gpuFunc.MSintegrateMult.MaxThreadsPerBlock 1 1];
else
    gpuFunc.active = false;
end

end