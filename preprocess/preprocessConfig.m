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

correlationActive = isfield(config.focalPointsL,'base_2');

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

% V shift
if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift'))
    [vShift_X,vShift_Y,vShift_Z] = ndgrid( ...
        config.focalPointsVshift.vShift, ...
        0, ...
        0);
    
    config.focalPointsVshift.vShiftVector = [vShift_X(:).'; vShift_Y(:).'; vShift_Z(:).'];
end

% Focal directions - L
if(isfield(config,'focalDirectionsL') && isfield(config.focalDirectionsL,'theta'))
    [T,P] = ndgrid( ...
        config.focalDirectionsL.theta, ...
        config.focalDirectionsL.phi);

    config.focalDirectionsL.vector = [ ...
        (sin(T(:)).').*(cos(P(:)).') ; ...
        (sin(T(:)).').*(sin(P(:)).') ; ...
        (cos(T(:)).') ]; 
elseif(isfield(config,'focalDirectionsL') && isfield(config.focalDirectionsL,'base'))
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
else
    config.focalDirectionsL.dim = config.focalPointsL.dim;
    config.focalDirectionsL.vector = 0 * config.focalPointsL.vector;
    config.focalDirectionsL.vector(end,:) = 1;
end

% Focal directions - V
if(isfield(config,'focalDirectionsV') && isfield(config.focalDirectionsV,'theta'))
    [T,P] = ndgrid( ...
        config.focalDirectionsV.theta, ...
        config.focalDirectionsV.phi);

    config.focalDirectionsV.vector = [ ...
        (sin(T(:)).').*(cos(P(:)).') ; ...
        (sin(T(:)).').*(sin(P(:)).') ; ...
        (cos(T(:)).') ]; 
elseif(isfield(config,'focalDirectionsV') && isfield(config.focalDirectionsV,'base'))
    if(config.focalDirectionsV.xyGrid)
        ldiryVec = sin(config.focalDirectionsV.base);
    else
        ldiryVec = 0;
    end

    [focalDirectionsV_X,focalDirectionsV_Y] = ndgrid( ...
        sin(config.focalDirectionsV.base), ...
        ldiryVec);

    config.focalDirectionsV.vector = [ ...
        focalDirectionsV_X(:).' ; ...
        focalDirectionsV_Y(:).' ; ...
        sqrt(1 - (focalDirectionsV_X(:).').^2 - (focalDirectionsV_Y(:).').^2)];
else
    config.focalDirectionsV.dim = config.focalPointsV.dim;
    config.focalDirectionsV.vector = 0 * config.focalPointsV.vector;
    config.focalDirectionsV.vector(end,:) = 1;
end

config.apertureVmf_l = movmfAperture( ...
    config.mask_varL , ...
    config.focalPointsL.vector , ...
    -1 , ...
    config.focalDirectionsL.vector , ...
    config.focalPointsL.dim , ...
    config.focalDirectionsL.dim);

if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift'))
    config.apertureVmf_v = movmfAperture( ...
        config.mask_varV , ...
        config.focalPointsV.vector , ...
        1 , ...
        config.focalDirectionsV.vector , ...
        config.focalPointsV.dim , ...
        config.focalDirectionsV.dim, ...
        config.focalPointsVshift.vShiftVector, ...
        config.focalPointsVshift.dim);
else
    config.apertureVmf_v = movmfAperture( ...
        config.mask_varV , ...
        config.focalPointsV.vector , ...
        1 , ...
        config.focalDirectionsV.vector , ...
        config.focalPointsV.dim , ...
        config.focalDirectionsV.dim);
end

%% Aperture - correlation
if(correlationActive)
    % Focal points - L

    if(config.focalPointsL.xyGrid)
        lyVec = config.focalPointsL.base_2;
    else
        lyVec = 0;
    end

    [focalPointsL_X,focalPointsL_Y,focalPointsL_Z] = ndgrid( ...
        config.focalPointsL.base_2, ...
        lyVec, ...
        config.focalPointsL.plain_2);
    config.focalPointsL.vector_2 = [focalPointsL_X(:).'; focalPointsL_Y(:).'; focalPointsL_Z(:).'];

    % Focal points - V
    if(config.focalPointsV.xyGrid)
        vyVec = config.focalPointsV.base_2;
    else
        vyVec = 0;
    end

    [focalPointsV_X,focalPointsV_Y,focalPointsV_Z] = ndgrid( ...
        config.focalPointsV.base_2, ...
        vyVec, ...
        config.focalPointsV.plain_2);
    config.focalPointsV.vector_2 = [focalPointsV_X(:).'; focalPointsV_Y(:).'; focalPointsV_Z(:).'];

    % V shift
    if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift_2'))
        [vShift_X,vShift_Y,vShift_Z] = ndgrid( ...
            config.focalPointsVshift.vShift_2, ...
            0, ...
            0);

        config.focalPointsVshift.vShiftVector_2 = [vShift_X(:).'; vShift_Y(:).'; vShift_Z(:).'];
    end

    % Focal directions - L
    if(isfield(config,'focalDirectionsL') && isfield(config.focalDirectionsL,'theta_2'))
        [T,P] = ndgrid( ...
            config.focalDirectionsL.theta_2, ...
            config.focalDirectionsL.phi_2);

        config.focalDirectionsL.vector_2 = [ ...
            (sin(T(:)).').*(cos(P(:)).') ; ...
            (sin(T(:)).').*(sin(P(:)).') ; ...
            (cos(T(:)).') ]; 
    elseif(isfield(config,'focalDirectionsL') && isfield(config.focalDirectionsL,'base_2'))
        if(config.focalDirectionsL.xyGrid)
            ldiryVec = sin(config.focalDirectionsL.base_2);
        else
            ldiryVec = 0;
        end

        [focalDirectionsL_X,focalDirectionsL_Y] = ndgrid( ...
            sin(config.focalDirectionsL.base_2), ...
            ldiryVec);

        config.focalDirectionsL.vector_2 = [ ...
            focalDirectionsL_X(:).' ; ...
            focalDirectionsL_Y(:).' ; ...
            sqrt(1 - (focalDirectionsL_X(:).').^2 - (focalDirectionsL_Y(:).').^2)];
    else
        config.focalDirectionsL.vector_2 = 0 * config.focalPointsL.vector_2;
        config.focalDirectionsL.vector_2(end,:) = 1;
    end

    % Focal directions - V
    if(isfield(config,'focalDirectionsV') && isfield(config.focalDirectionsV,'theta_2'))
        [T,P] = ndgrid( ...
            config.focalDirectionsV.theta_2, ...
            config.focalDirectionsV.phi_2);

        config.focalDirectionsV.vector_2 = [ ...
            (sin(T(:)).').*(cos(P(:)).') ; ...
            (sin(T(:)).').*(sin(P(:)).') ; ...
            (cos(T(:)).') ]; 
    elseif(isfield(config,'focalDirectionsV') && isfield(config.focalDirectionsV,'base_2'))
        if(config.focalDirectionsV.xyGrid)
            ldiryVec = sin(config.focalDirectionsV.base_2);
        else
            ldiryVec = 0;
        end

        [focalDirectionsV_X,focalDirectionsV_Y] = ndgrid( ...
            sin(config.focalDirectionsV.base_2), ...
            ldiryVec);

        config.focalDirectionsV.vector_2 = [ ...
            focalDirectionsV_X(:).' ; ...
            focalDirectionsV_Y(:).' ; ...
            sqrt(1 - (focalDirectionsV_X(:).').^2 - (focalDirectionsV_Y(:).').^2)];
    else
        config.focalDirectionsV.vector_2 = 0 * config.focalPointsV.vector_2;
        config.focalDirectionsV.vector_2(end,:) = 1;
    end

    apertureVmf_l_2 = movmfAperture( ...
        config.mask_varL , ...
        config.focalPointsL.vector_2 , ...
        -1 , ...
        config.focalDirectionsL.vector_2 , ...
        config.focalPointsL.dim , ...
        config.focalDirectionsL.dim);

    if(isfield(config,'focalPointsVshift') && isfield(config.focalPointsVshift,'vShift_2'))
        apertureVmf_v_2 = movmfAperture( ...
            config.mask_varV , ...
            config.focalPointsV.vector_2 , ...
            1 , ...
            config.focalDirectionsV.vector_2 , ...
            config.focalPointsV.dim , ...
            config.focalDirectionsV.dim, ...
            config.focalPointsVshift.vShiftVector_2, ...
            config.focalPointsVshift.dim);
    else
        apertureVmf_v_2 = movmfAperture( ...
            config.mask_varV , ...
            config.focalPointsV.vector_2 , ...
            1 , ...
            config.focalDirectionsV.vector_2 , ...
            config.focalPointsV.dim , ...
            config.focalDirectionsV.dim);
    end

    config.apertureVmf_l = movmfUnite(config.apertureVmf_l,apertureVmf_l_2);
    config.apertureVmf_v = movmfUnite(config.apertureVmf_v,apertureVmf_v_2);
    
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
pxItersNum = 1e3;
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
