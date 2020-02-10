function config = preprocessConfig_onGrid(config)
%% Sample
config.box_min = [-config.boxAxial/2;-config.boxAxial/2;-config.boxDepth/2];
config.box_max = [config.boxAxial/2;config.boxAxial/2;config.boxDepth/2];

if(config.mcGpuL)
    config.box_min = gpuArray(config.box_min);
    config.box_max = gpuArray(config.box_max);
end

%% Aperture
% Focal points
[focalPoints_X,focalPoints_Y,focalPoints_Z] = ndgrid( ...
    config.focalPoints_base, ...
    config.focalPoints_base, ...
    config.focalPoints_plain);

config.focalPointsL = [focalPoints_X(:).'; ...
            focalPoints_Y(:).'; ...
            focalPoints_Z(:).' ];
        
config.focalPointsV = config.focalPointsL;
        
% Focal directions
[focalLDirections_X,focalLDirections_Y] = ndgrid( ...
    config.focalLDirections, ...
    config.focalLDirections);

config.focalDirectionsL = [ ...
        sin(focalLDirections_X(:).'); ...
        sin(focalLDirections_Y(:).') ; ...
        sqrt(1 - sin(focalLDirections_X(:).').^2 - sin(focalLDirections_Y(:).').^2 ) ];

[focalVDirections_X,focalVDirections_Y] = ndgrid( ...
    config.focalVDirections, ...
    config.focalVDirections);

config.focalDirectionsV = [ ...
        sin(focalVDirections_X(:).'); ...
        sin(focalVDirections_Y(:).') ; ...
        sqrt(1 - sin(focalVDirections_X(:).').^2 - sin(focalVDirections_Y(:).').^2 ) ];

if(config.mcGpuL)
    config.focalPointsL = gpuArray(config.focalPointsL);
    config.focalDirectionsL = gpuArray(config.focalDirectionsL);
end

if(config.mcGpuV)
    config.focalPointsV = gpuArray(config.focalPointsV);
    config.focalDirectionsV = gpuArray(config.focalDirectionsV);
end

config.mask_varL = config.mask_var;
config.mask_varV = config.mask_var;
config.apertureVmf_l = movmfAperture(config.mask_var,config.focalPointsL,-1,config.focalDirectionsL,2,3);
config.apertureVmf_v = movmfAperture(config.mask_var,config.focalPointsV,1,config.focalDirectionsV,2,4);

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
    config.smpPreprocess = preprocess_smpVmfBeamSum(config,1e2);
end

px = zeros(1,pxItersNum);

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

end