function [config] = preprocess_near_field(config)
%% Simluation
config.simulation.gpuNum = uint32(config.simulation.gpuNum);
config.simulation.iterations = uint32(config.simulation.iterations);

if(~isfield(config.simulation,'precision'))
    config.simulation.precision = 'double';
end

%% Sample
config.medium.box_min = [-config.medium.boxAxial/2;-config.medium.boxAxial/2;-config.medium.boxDepth/2];
config.medium.box_max = [config.medium.boxAxial/2;config.medium.boxAxial/2;config.medium.boxDepth/2];

if(isfield(config,'boxShift'))
    config.medium.box_min = config.medium.box_min + config.medium.boxShift;
    config.medium.box_max = config.medium.box_max + config.medium.boxShift;
end

config.medium.sigt = 1/config.medium.MFP;
correlationActive = isfield(config.nf.focalPointsL,'x_2');

%% Parameters
currDim = 1;

paramsList = fieldnames(config.nf.parameters);

for paramNum = 1:1:numel(paramsList)
    param = config.nf.parameters.(paramsList{paramNum});
    param = param(:);
    permuteVec = 1:1:currDim;
    permuteVec(1) = currDim; permuteVec(end) = 1;
    if(numel(permuteVec) > 1)
        param = permute(param,permuteVec);
    end
    
    config.nf.parameters.(paramsList{paramNum}) = param;
    currDim = currDim + 1;
end

%% Aperture

if(~isfield(config.aperture,'is_normalized'))
    config.aperture.is_normalized = true;
end

config.aperture.kappa_l = 1/(config.aperture.mask_varL^2);
config.aperture.kappa_v = 1/(config.aperture.mask_varV^2);

if(config.aperture.is_normalized)
    config.aperture.c_l = ((3/2-1)*log(config.aperture.kappa_l) - (3/2)*log(2*pi) - logbesseli(3,config.aperture.kappa_l));
    config.aperture.c_v = ((3/2-1)*log(config.aperture.kappa_v) - (3/2)*log(2*pi) - logbesseli(3,config.aperture.kappa_v));
else
    config.aperture.c_l = -config.aperture.kappa_l;
    config.aperture.c_v = -config.aperture.kappa_v;
end


%% Evaluate

if(~isfield(config.medium,'wavelenght'))
    config.nf.wavenumber_st.k = eval(['config.nf.wavenumber.k',parseFunc(config.nf.wavenumber.k,'config.nf.parameters.')]);
else
    config.nf.wavenumber_st.k = config.medium.wavelenght;
end

config.nf.l_st.P1 = eval(['config.nf.focalPointsL.x',parseFunc(config.nf.focalPointsL.x)]);
config.nf.l_st.P2 = eval(['config.nf.focalPointsL.y',parseFunc(config.nf.focalPointsL.y)]);
config.nf.l_st.P3 = eval(['config.nf.focalPointsL.z',parseFunc(config.nf.focalPointsL.z)]);

config.nf.v_st.P1 = eval(['config.nf.focalPointsV.x',parseFunc(config.nf.focalPointsV.x)]);
config.nf.v_st.P2 = eval(['config.nf.focalPointsV.y',parseFunc(config.nf.focalPointsV.y)]);
config.nf.v_st.P3 = eval(['config.nf.focalPointsV.z',parseFunc(config.nf.focalPointsV.z)]);

config.nf.l_st.D1 = eval(['config.nf.focalDirectionsL.x',parseFunc(config.nf.focalDirectionsL.x)]);
config.nf.l_st.D2 = eval(['config.nf.focalDirectionsL.y',parseFunc(config.nf.focalDirectionsL.y)]);
config.nf.l_st.D3 = eval(['config.nf.focalDirectionsL.z',parseFunc(config.nf.focalDirectionsL.z)]);

config.nf.v_st.D1 = eval(['config.nf.focalDirectionsV.x',parseFunc(config.nf.focalDirectionsV.x)]);
config.nf.v_st.D2 = eval(['config.nf.focalDirectionsV.y',parseFunc(config.nf.focalDirectionsV.y)]);
config.nf.v_st.D3 = eval(['config.nf.focalDirectionsV.z',parseFunc(config.nf.focalDirectionsV.z)]);

if(correlationActive)
    if(~isfield(config.medium,'wavelenght'))
        config.nf.wavenumber_st.k_2 = eval(['config.nf.wavenumber.k_2',parseFunc(config.nf.wavenumber.k_2,'config.nf.parameters.')]);
    else
        config.nf.wavenumber_st.k_2 = config.medium.wavelenght;
    end

    config.nf.l_st.P1_2 = eval(['config.nf.focalPointsL.x_2',parseFunc(config.nf.focalPointsL.x_2)]);
    config.nf.l_st.P2_2 = eval(['config.nf.focalPointsL.y_2',parseFunc(config.nf.focalPointsL.y_2)]);
    config.nf.l_st.P3_2 = eval(['config.nf.focalPointsL.z_2',parseFunc(config.nf.focalPointsL.z_2)]);

    config.nf.v_st.P1_2 = eval(['config.nf.focalPointsV.x_2',parseFunc(config.nf.focalPointsV.x_2)]);
    config.nf.v_st.P2_2 = eval(['config.nf.focalPointsV.y_2',parseFunc(config.nf.focalPointsV.y_2)]);
    config.nf.v_st.P3_2 = eval(['config.nf.focalPointsV.z_2',parseFunc(config.nf.focalPointsV.z_2)]);

    config.nf.l_st.D1_2 = eval(['config.nf.focalDirectionsL.x_2',parseFunc(config.nf.focalDirectionsL.x_2)]);
    config.nf.l_st.D2_2 = eval(['config.nf.focalDirectionsL.y_2',parseFunc(config.nf.focalDirectionsL.y_2)]);
    config.nf.l_st.D3_2 = eval(['config.nf.focalDirectionsL.z_2',parseFunc(config.nf.focalDirectionsL.z_2)]);

    config.nf.v_st.D1_2 = eval(['config.nf.focalDirectionsV.x_2',parseFunc(config.nf.focalDirectionsV.x_2)]);
    config.nf.v_st.D2_2 = eval(['config.nf.focalDirectionsV.y_2',parseFunc(config.nf.focalDirectionsV.y_2)]);
    config.nf.v_st.D3_2 = eval(['config.nf.focalDirectionsV.z_2',parseFunc(config.nf.focalDirectionsV.z_2)]);
end

%% Destroy all function handlers
fn = fieldnames(config.nf.wavenumber);
for k=1:numel(fn)
    config.nf.wavenumber.(fn{k}) = func2str(config.nf.wavenumber.(fn{k}));
end

fn = fieldnames(config.nf.focalPointsL);
for k=1:numel(fn)
    config.nf.focalPointsL.(fn{k}) = func2str(config.nf.focalPointsL.(fn{k}));
end

fn = fieldnames(config.nf.focalPointsV);
for k=1:numel(fn)
    config.nf.focalPointsV.(fn{k}) = func2str(config.nf.focalPointsV.(fn{k}));
end

fn = fieldnames(config.nf.focalDirectionsL);
for k=1:numel(fn)
    config.nf.focalDirectionsL.(fn{k}) = func2str(config.nf.focalDirectionsL.(fn{k}));
end

fn = fieldnames(config.nf.focalDirectionsV);
for k=1:numel(fn)
    config.nf.focalDirectionsV.(fn{k}) = func2str(config.nf.focalDirectionsV.(fn{k}));
end

%% Scattering function
if(config.scatter.type == 3)
    if(~isfield(config.scatter,'forwardWeight'))
        config.scatter.forwardWeight = 1;
    end
end

%% Scattering function - entry
entry.sctType = config.scatter.type;
entry.dim = config.medium.dim;
entry.vmf_k = config.movmf.vmf_k;
entry.vmf_iterations = config.movmf.vmf_iterations;

if(config.scatter.type == 3)
    entry.g = config.scatter.g;
    entry.forwardWeight = config.scatter.forwardWeight;
    entry.vmf_samples = config.movmf.vmf_samples;
    vmf_samples = entry.vmf_samples;
end

if(config.scatter.type == 2)
    entry.pdf = config.scatter.pdf;
    vmf_samples = 0;
end

%% Scattering function - throughput
entryHash = DataHash(entry);

preprocessPath = which('vmfCache.mat');
if(isempty(preprocessPath))
    if(config.scatter.type == 3)
        ampfunc.g = config.scatter.g;
        ampfunc.forwardWeight = config.scatter.forwardWeight;
    end

    if(config.scatter.type == 2)
        theta = linspace(0,pi,numel(config.scatter.pdf));
        f = config.scatter.pdf;
        
        ampfunc.cs = mean((abs(f).^2).*sin(theta))*2*pi^2;
        ampfunc.samplePdf = (abs(f(:).').^2) .* sin(theta) ./ ...
            sum((abs(f(:).').^2) .* sin(theta));
        ampfunc.sampleCdf = cumsum(config.ampfunc.samplePdf);
        
        ampfunc.evalPdf = (abs(f(:).').^2) ./ config.ampfunc.cs;
        ampfunc.evalAmp = config.ampfunc.evalPdf .^ 0.5;
        
        ampfunc.sampleIcdf = invCDFvec(config.ampfunc.sampleCdf);
    end
    
    config.movmf = movmfBuild(config.scatter.type, ampfunc, config.medium.dim, ...
        config.movmf.vmf_k, vmf_samples, config.movmf.vmf_iterations, true);
    
    vmfCache(1).entryHash = entryHash;
    vmfCache(1).entry = entry;
    vmfCache(1).ampfunc = ampfunc;
    vmfCache(1).movmf = config.movmf;
    
    config.movmf.cacheIdx = 1;
    
    currentFilePath = which('preprocess_near_field.m');
    currentFileFolder = regexp(currentFilePath,'**[\/\\]preprocess_near_field.m');
    
    save([currentFilePath(1:currentFileFolder),'vmfCache.mat'],'vmfCache');
else
    load(preprocessPath,'vmfCache')
    
    isInCache = false;
    for idxNum = 1:1:numel(vmfCache)
        if(strcmp(entryHash,vmfCache(idxNum).entryHash))
            config.movmf = vmfCache(idxNum).movmf;
            config.movmf.cacheIdx = idxNum;
            isInCache = true;
            break;
        end
    end
    
    if(~isInCache)
        lastElem = numel(vmfCache);
        
        if(config.scatter.type == 3)
            ampfunc.g = config.scatter.g;
            ampfunc.forwardWeight = config.scatter.forwardWeight;
        end

        if(config.scatter.type == 2)
            theta = linspace(0,pi,numel(config.scatter.pdf));
            f = config.scatter.pdf;

            ampfunc.cs = mean((abs(f).^2).*sin(theta))*2*pi^2;
            ampfunc.samplePdf = (abs(f(:).').^2) .* sin(theta) ./ ...
                sum((abs(f(:).').^2) .* sin(theta));
            ampfunc.sampleCdf = cumsum(config.ampfunc.samplePdf);

            ampfunc.evalPdf = (abs(f(:).').^2) ./ config.ampfunc.cs;
            ampfunc.evalAmp = config.ampfunc.evalPdf .^ 0.5;

            ampfunc.sampleIcdf = invCDFvec(config.ampfunc.sampleCdf);
        end
        
        config.movmf = movmfBuild(config.scatter.type, ampfunc, config.medium.dim, ...
            config.movmf.vmf_k, vmf_samples, config.movmf.vmf_iterations, true);
        
        vmfCache(lastElem+1).entryHash = entryHash;
        vmfCache(lastElem+1).entry = entry;
        vmfCache(lastElem+1).ampfunc = ampfunc;
        vmfCache(lastElem+1).movmf = config.movmf;
        config.movmf.cacheIdx = lastElem + 1;
        
        save(preprocessPath,'vmfCache');
    end
end

%% Importance sampling
config.sample.position_type = uint32(config.sample.position_type);
config.sample.direction_type = uint32(config.sample.direction_type);

if(~isfield(config.sample,'same_beam'))
    config.sample.same_beam = true;
end

config.sample.same_beam = uint32(config.sample.same_beam);

if(~isfield(config.sample,'z0_samples'))
    if(correlationActive)
        config.sample.z0_samples = unique([config.nf.l_st.P3(:);config.nf.l_st.P3_2(:)],'sorted');
    else
        config.sample.z0_samples = unique(config.nf.l_st.P3(:),'sorted');
    end
end

config.sample.z0_sample_num = uint32(numel(config.sample.z0_samples));

if(~isfield(config.sample,'z_sample_num'))
    config.sample.z_sample_num = uint32(10);
end

config.sample.z_sample_num = uint32(config.sample.z_sample_num);

if(isfield(config.sample,'min_px0') || isfield(config.sample,'min_pw0'))
    config.sample.is_min_known = uint32(true);
    
    if(isfield(config.sample,'test_rounds'))
        config.sample = rmfield(config.sample,'test_rounds');
    end
    
    if(isfield(config.sample,'min_probability_percent'))
        config.sample = rmfield(config.sample,'min_probability_percent');
    end
    
    if(~isfield(config.sample,'min_px0'))
        config.sample.min_px0 = 0.0;
    end
    
    if(~isfield(config.sample,'min_pw0'))
        config.sample.min_pw0 = 0.0;
    end
else
    config.sample.is_min_known = uint32(false);
    
    if(isfield(config.sample,'test_rounds'))
        config.sample.test_rounds = uint32(config.sample.test_rounds);
    else
        config.sample.test_rounds = uint32(10);
    end
    
    if(~isfield(config.sample,'min_probability_percent'))
        config.sample.min_probability_percent = 0.01;
    end
end

%% Parpool
if(numel(config.simulation.gpuNum) > 1)
    if(isempty(gcp('nocreate')))
        parpool(numel(config.gpuNum));
    end
end

%% precision

config.scatter.type = uint32(config.scatter.type);
config.medium.dim = uint32(config.medium.dim);

if(strcmp(config.simulation.precision,'single'))
    config.medium.box_min = single(config.medium.box_min);
    config.medium.box_max = single(config.medium.box_max);
    config.medium.sigt = single(config.medium.sigt);
    config.movmf.alpha = single(config.movmf.alpha);
    config.movmf.mu3 = single(config.movmf.mu3);
    config.movmf.c = single(config.movmf.c);
    
    config.aperture.kappa_l = single(config.aperture.kappa_l);
    config.aperture.c_l = single(config.aperture.c_l);
    config.aperture.kappa_v = single(config.aperture.kappa_v);
    config.aperture.c_v = single(config.aperture.c_v);
    
    config.nf.wavenumber_st.k = single(config.nf.wavenumber_st.k);
    config.nf.l_st.P1 = single(config.nf.l_st.P1);
    config.nf.l_st.P2 = single(config.nf.l_st.P2);
    config.nf.l_st.P3 = single(config.nf.l_st.P3);
    config.nf.l_st.D1 = single(config.nf.l_st.D1);
    config.nf.l_st.D2 = single(config.nf.l_st.D2);
    config.nf.l_st.D3 = single(config.nf.l_st.D3);
    config.nf.v_st.P1 = single(config.nf.v_st.P1);
    config.nf.v_st.P2 = single(config.nf.v_st.P2);
    config.nf.v_st.P3 = single(config.nf.v_st.P3);
    config.nf.v_st.D1 = single(config.nf.v_st.D1);
    config.nf.v_st.D2 = single(config.nf.v_st.D2);
    config.nf.v_st.D3 = single(config.nf.v_st.D3);
    
    config.sample.z0_samples = single(config.sample.z0_samples);
    
    if(isfield(config.sample,'min_px0'))
        config.sample.min_px0 = single(config.sample.min_px0);
    end
    
    if(isfield(config.sample,'min_pw0'))
        config.sample.min_pw0 = single(config.sample.min_pw0);
    end
    
    if(isfield(config.sample,'min_probability_percent'))
        config.sample.min_probability_percent = single(config.sample.min_probability_percent);
    end
    
    if(config.scatter.type == 3)
        config.scatter.g = single(config.scatter.g);
    end
    
    if(correlationActive)
        config.nf.wavenumber_st.k_2 = single(config.nf.wavenumber_st.k_2);
        config.nf.l_st.P1_2 = single(config.nf.l_st.P1_2);
        config.nf.l_st.P2_2 = single(config.nf.l_st.P2_2);
        config.nf.l_st.P3_2 = single(config.nf.l_st.P3_2);
        config.nf.l_st.D1_2 = single(config.nf.l_st.D1_2);
        config.nf.l_st.D2_2 = single(config.nf.l_st.D2_2);
        config.nf.l_st.D3_2 = single(config.nf.l_st.D3_2);
        config.nf.v_st.P1_2 = single(config.nf.v_st.P1_2);
        config.nf.v_st.P2_2 = single(config.nf.v_st.P2_2);
        config.nf.v_st.P3_2 = single(config.nf.v_st.P3_2);
        config.nf.v_st.D1_2 = single(config.nf.v_st.D1_2);
        config.nf.v_st.D2_2 = single(config.nf.v_st.D2_2);
        config.nf.v_st.D3_2 = single(config.nf.v_st.D3_2);
    end
end

end
