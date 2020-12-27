function [config] = preprocess_refocus(config)
%% Simluation
if(config.refocus.sample_random)
    config.simulation.iterations = config.simulation.iterations * 1024;
end

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

%% Refocus
if(~isfield(config.refocus,'sample_random'))
    config.refocus.sample_random = false;
end

if(~isfield(config.refocus,'binary_aperture'))
    config.refocus.binary_aperture = false;
end

if(~isfield(config.refocus,'sample_forward'))
    config.refocus.sample_forward = true;
end

if(~isfield(config.refocus,'sample_backward'))
    config.refocus.sample_backward = true;
end

if(~isfield(config.refocus,'max_xy_value'))
    config.refocus.max_xy_value = 1;
end

if(~isfield(config.refocus,'bias_attenuation'))
    config.refocus.bias_attenuation = false;
end

if(~isfield(config.refocus,'tabulated_dldv'))
    config.refocus.tabulated_dldv = 1/config.medium.boxAxial;
end


if(~isfield(config.refocus,'random_directions_number'))
    config.refocus.random_directions_number = 64;
end

config.refocus.sample_random = uint32(config.refocus.sample_random);
config.refocus.binary_aperture = uint32(config.refocus.binary_aperture);
config.refocus.sample_forward = uint32(config.refocus.sample_forward);
config.refocus.sample_backward = uint32(config.refocus.sample_backward);
config.refocus.bias_attenuation = uint32(config.refocus.bias_attenuation);
config.refocus.random_directions_number = uint32(config.refocus.random_directions_number);

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
