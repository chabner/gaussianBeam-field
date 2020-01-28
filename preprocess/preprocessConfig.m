function config = preprocessConfig(config)
%% Sample
if(config.dimNum == 2)
    config.box_min = [-config.boxAxial/2;-config.boxDepth/2];
    config.box_max = [config.boxAxial/2;config.boxDepth/2];
else
    config.box_min = [-config.boxAxial/2;-config.boxAxial/2;-config.boxDepth/2];
    config.box_max = [config.boxAxial/2;config.boxAxial/2;config.boxDepth/2];
end

%% Aperture
% Focal points
config.focalPointsL = [config.focalPointsL_base; ...
            0 * config.focalPointsL_base; ...
            config.focalPointsL_plain * (config.focalPointsL_base.^0) ];

[focalPointsV_X,focalPointsV_Y] = ndgrid(config.focalPointsV_base,config.focalPointsV_base);

config.focalPointsV = [focalPointsV_X(:).'; ...
            focalPointsV_Y(:).'; ...
            config.focalPointsV_plain * ones(1,numel(focalPointsV_X)) ];
        
% Focal directions
if(isfield(config,'focalLDirections'))
    config.focalDirectionsL = [ ...
        sin(config.focalLDirections); ...
        0 * config.focalLDirections ; ...
        cos(config.focalLDirections) ];
else
    config.focalDirectionsL = [ ...
        0 * config.focalLDirections ; ...
        0 * config.focalLDirections ; ...
        config.focalLDirections .^ 0 ];
end

if(isfield(config,'focalVDirections'))
    config.focalDirectionsV = [ ...
        sin(config.focalVDirections); ...
        0 * config.focalVDirections ; ...
        cos(config.focalVDirections) ];
else
    config.focalDirectionsV = [ ...
        0 * config.focalVDirections ; ...
        0 * config.focalVDirections ; ...
        config.focalVDirections .^ 0 ];
end

%% Scattering function
if(config.sctType == 3)
    config.sct_type = 3;
    config.ampfunc.g = config.g;
    config.ampfunc.forwardWeight = config.forwardWeight;
    config.ampfunc0 = inf;
end

%% Scattering function
if(~isfile(['preprocess',filesep,'vmfCache.mat']))
    config.movmf = movmfBuild(config.sctType, config.ampfunc, config.dimNum, config.vmf_k, config.vmf_samples, config.vmf_iterations, config.useGpu);
    
    vmfCache(1).g = config.g;
    vmfCache(1).forwardWeight = config.forwardWeight;
    vmfCache(1).vmf_k = config.vmf_k;
    vmfCache(1).vmf_iterations = config.vmf_iterations;
    vmfCache(1).vmf_samples = config.vmf_samples;
    vmfCache(1).movmf = config.movmf;
    
    save(['preprocess',filesep,'vmfCache.mat'],'vmfCache');
else
    load(['preprocess',filesep,'vmfCache.mat'],'vmfCache')
    
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
    
        save(['preprocess',filesep,'vmfCache.mat'],'vmfCache');
    end
end

%% Importance sampling
% Sample 100 different samples, and take the median of px, in order to
% estimate bad samples
pxItersNum = 100;
minTorr = 1e-2;

box_w = config.box_max-config.box_min;
V=prod(box_w);

if(mod(config.sampleFlag,10) == 3)
    config.smpPreprocess = preprocess_smpVmfBeam(config,lightNum,1e4);
end

if(mod(config.sampleFlag,10) == 4)
    config.smpPreprocess = preprocess_smpVmfBeamSum(config,1e4);
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
        apertureVmf_l = movmfAperture(config.mask_varL,config.focalPointsL,-1,config.focalDirectionsL);
        [~,px(iterNum)] = smpVmfBeamSum(apertureVmf_l,config.smpPreprocess,config.box_min,config.box_max);
        px(iterNum) = px(iterNum) * V;
    end
end

config.smpPreprocess(1).pxMin = median(px) * minTorr;

end