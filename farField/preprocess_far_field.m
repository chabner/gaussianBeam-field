function [config] = preprocess_far_field(config)
%% Simulation
config.simulation.gpuNum = uint32(config.simulation.gpuNum);
config.simulation.iterations = uint32(config.simulation.iterations);

if(~isfield(config.simulation,'cbs'))
    config.simulation.cbs = false;
end

config.simulation.cbs = uint32(config.simulation.cbs);

if(~isfield(config.simulation,'precision'))
    config.simulation.precision = 'double';
end

%% Sample
if(config.medium.dim == 3)
    config.medium.box_min = [-config.medium.boxAxial/2;-config.medium.boxAxial/2;-config.medium.boxDepth/2];
    config.medium.box_max = [config.medium.boxAxial/2;config.medium.boxAxial/2;config.medium.boxDepth/2];
end

if(config.medium.dim == 2)
    config.medium.box_min = [-config.medium.boxAxial/2;-config.medium.boxDepth/2];
    config.medium.box_max = [config.medium.boxAxial/2;config.medium.boxDepth/2];
end

if(isfield(config.medium,'boxShift'))
    config.medium.box_min = config.medium.box_min + config.medium.boxShift;
    config.medium.box_max = config.medium.box_max + config.medium.boxShift;
end

config.medium.sigt = 1/config.medium.MFP;

noAttDir = ~isfield(config.ff,'dir_v');
ff_correlationActive = isfield(config.ff.l,'x_2');

%% Parameters
currDim = 1;

% ff parameters
paramsList = fieldnames(config.ff.parameters);
ff_dims = numel(paramsList);

for paramNum = 1:1:ff_dims
    param = config.ff.parameters.(paramsList{paramNum});
    param = param(:);
    permuteVec = 1:1:currDim;
    permuteVec(1) = currDim; permuteVec(end) = 1;
    if(numel(permuteVec) > 1)
        param = permute(param,permuteVec);
    end

    config.ff.parameters.(paramsList{paramNum}) = param;
    currDim = currDim + 1;
end

%% Evaluate
if(~isfield(config.medium,'wavelenght'))
    config.ff.wavenumber_st.k = eval(['config.ff.wavenumber.k',parseFunc(config.ff.wavenumber.k,'config.ff.parameters.')]);
else
    config.ff.wavenumber_st.k = config.medium.wavelenght;
end


config.ff.v_st.X = eval(['config.ff.v.x',parseFunc(config.ff.v.x,'config.ff.parameters.')]);
config.ff.v_st.Y = eval(['config.ff.v.y',parseFunc(config.ff.v.y,'config.ff.parameters.')]);

if(config.medium.dim == 3)
    config.ff.v_st.Z = eval(['config.ff.v.z',parseFunc(config.ff.v.z,'config.ff.parameters.')]);
else
    config.ff.v_st.Z = 0;
end

config.ff.l_st.X = eval(['config.ff.l.x',parseFunc(config.ff.l.x,'config.ff.parameters.')]);
config.ff.l_st.Y = eval(['config.ff.l.y',parseFunc(config.ff.l.y,'config.ff.parameters.')]);

if(config.medium.dim == 3)
    config.ff.l_st.Z = eval(['config.ff.l.z',parseFunc(config.ff.l.z,'config.ff.parameters.')]);
else
    config.ff.l_st.Z = 0;
end

if(~noAttDir)
    config.ff.v_st.DIR_X = eval(['config.ff.dir_v.x',parseFunc(config.ff.dir_v.x,'config.ff.parameters.')]);
    config.ff.v_st.DIR_Y = eval(['config.ff.dir_v.y',parseFunc(config.ff.dir_v.y,'config.ff.parameters.')]);
    
    if(config.medium.dim == 3)
        config.ff.v_st.DIR_Z = eval(['config.ff.dir_v.z',parseFunc(config.ff.dir_v.z,'config.ff.parameters.')]);
    else
        config.ff.v_st.DIR_Z = 0;
    end

    config.ff.l_st.DIR_X = eval(['config.ff.dir_l.x',parseFunc(config.ff.dir_l.x,'config.ff.parameters.')]);
    config.ff.l_st.DIR_Y = eval(['config.ff.dir_l.y',parseFunc(config.ff.dir_l.y,'config.ff.parameters.')]);
    
    if(config.medium.dim == 3)
        config.ff.l_st.DIR_Z = eval(['config.ff.dir_l.z',parseFunc(config.ff.dir_l.z,'config.ff.parameters.')]);
    else
        config.ff.l_st.DIR_Z = 0;
    end
else
    config.ff.v_st.DIR_X = config.ff.v_st.X;
    config.ff.v_st.DIR_Y = config.ff.v_st.Y;
    
    if(config.medium.dim == 3)
        config.ff.v_st.DIR_Z = config.ff.v_st.Z;
    else
        config.ff.v_st.DIR_Z = 0;
    end

    config.ff.l_st.DIR_X = config.ff.l_st.X;
    config.ff.l_st.DIR_Y = config.ff.l_st.Y;
    
    if(config.medium.dim == 3)
        config.ff.l_st.DIR_Z = config.ff.l_st.Z;
    else
        config.ff.l_st.DIR_Z = 0;
    end
end

if(ff_correlationActive)
    if(~isfield(config.medium,'wavelenght'))
        config.ff.wavenumber_st.k_2 = eval(['config.ff.wavenumber.k_2',parseFunc(config.ff.wavenumber.k_2,'config.ff.parameters.')]);
    else
        config.ff.wavenumber_st.k_2 = config.medium.wavelenght;
    end

    config.ff.v_st.X_2 = eval(['config.ff.v.x_2',parseFunc(config.ff.v.x_2,'config.ff.parameters.')]);
    config.ff.v_st.Y_2 = eval(['config.ff.v.y_2',parseFunc(config.ff.v.y_2,'config.ff.parameters.')]);
    
    if(config.medium.dim == 3)
        config.ff.v_st.Z_2 = eval(['config.ff.v.z_2',parseFunc(config.ff.v.z_2,'config.ff.parameters.')]);
    else
        config.ff.v_st.Z_2 = 0;
    end

    config.ff.l_st.X_2 = eval(['config.ff.l.x_2',parseFunc(config.ff.l.x_2,'config.ff.parameters.')]);
    config.ff.l_st.Y_2 = eval(['config.ff.l.y_2',parseFunc(config.ff.l.y_2,'config.ff.parameters.')]);
    
    if(config.medium.dim == 3)
        config.ff.l_st.Z_2 = eval(['config.ff.l.z_2',parseFunc(config.ff.l.z_2,'config.ff.parameters.')]);
    else
        config.ff.l_st.Z_2 = 0;
    end

    if(~noAttDir)
        config.ff.v_st.DIR_X_2 = eval(['config.ff.dir_v.x_2',parseFunc(config.ff.dir_v.x_2,'config.ff.parameters.')]);
        config.ff.v_st.DIR_Y_2 = eval(['config.ff.dir_v.y_2',parseFunc(config.ff.dir_v.y_2,'config.ff.parameters.')]);
        
        if(config.medium.dim == 3)
            config.ff.v_st.DIR_Z_2 = eval(['config.ff.dir_v.z_2',parseFunc(config.ff.dir_v.z_2,'config.ff.parameters.')]);
        else
            config.ff.v_st.DIR_Z_2 = 0;
        end

        config.ff.l_st.DIR_X_2 = eval(['config.ff.dir_l.x_2',parseFunc(config.ff.dir_l.x_2,'config.ff.parameters.')]);
        config.ff.l_st.DIR_Y_2 = eval(['config.ff.dir_l.y_2',parseFunc(config.ff.dir_l.y_2,'config.ff.parameters.')]);
        
        if(config.medium.dim == 3)
            config.ff.l_st.DIR_Z_2 = eval(['config.ff.dir_l.z_2',parseFunc(config.ff.dir_l.z_2,'config.ff.parameters.')]);
        else
            config.ff.l_st.DIR_Z_2 = 0;
        end
    else
        config.ff.v_st.DIR_X_2 = config.ff.v_st.X_2;
        config.ff.v_st.DIR_Y_2 = config.ff.v_st.Y_2;
        
        if(config.medium.dim == 3)
            config.ff.v_st.DIR_Z_2 = config.ff.v_st.Z_2;
        else
            config.ff.v_st.DIR_Z_2 = 0;
        end

        config.ff.l_st.DIR_X_2 = config.ff.l_st.X_2;
        config.ff.l_st.DIR_Y_2 = config.ff.l_st.Y_2;
        
        if(config.medium.dim == 3)
            config.ff.l_st.DIR_Z_2 = config.ff.l_st.Z_2;
        else
            config.ff.l_st.DIR_Z_2 = 0;
        end
    end
end

%% Destroy all function handlers
fn = fieldnames(config.ff.wavenumber);
for k=1:numel(fn)
    config.ff.wavenumber.(fn{k}) = func2str(config.ff.wavenumber.(fn{k}));
end

fn = fieldnames(config.ff.v);
for k=1:numel(fn)
    config.ff.v.(fn{k}) = func2str(config.ff.v.(fn{k}));
end

fn = fieldnames(config.ff.l);
for k=1:numel(fn)
    config.ff.l.(fn{k}) = func2str(config.ff.l.(fn{k}));
end

%% Sampling
config.sample.position_type = uint32(config.sample.position_type);
config.sample.direction_type = uint32(config.sample.direction_type);

if(~isfield(config.sample,'mean_l'))
    config.sample.mean_l = true;
end

config.sample.mean_l = uint32(config.sample.mean_l);

if(config.sample.direction_type == 2)
    if(isfield(config.sample,'f0'))
        config.sample.f0 = complex(config.sample.f0);
    end
end

%% Precision
if(strcmp(config.simulation.precision,'single'))
    config.medium.box_min = single(config.medium.box_min);
    config.medium.box_max = single(config.medium.box_max);
    config.medium.sigt = single(config.medium.sigt);
    config.ff.wavenumber_st.k = single(config.ff.wavenumber_st.k);
    config.ff.l_st.X = single(config.ff.l_st.X);
    config.ff.l_st.Y = single(config.ff.l_st.Y);
    config.ff.l_st.Z = single(config.ff.l_st.Z);
    config.ff.l_st.DIR_X = single(config.ff.l_st.DIR_X);
    config.ff.l_st.DIR_Y = single(config.ff.l_st.DIR_Y);
    config.ff.l_st.DIR_Z = single(config.ff.l_st.DIR_Z);
    config.ff.v_st.X = single(config.ff.v_st.X);
    config.ff.v_st.Y = single(config.ff.v_st.Y);
    config.ff.v_st.Z = single(config.ff.v_st.Z);
    config.ff.v_st.DIR_X = single(config.ff.v_st.DIR_X);
    config.ff.v_st.DIR_Y = single(config.ff.v_st.DIR_Y);
    config.ff.v_st.DIR_Z = single(config.ff.v_st.DIR_Z);
    
    if(ff_correlationActive)
        config.ff.wavenumber_st.k_2 = single(config.ff.wavenumber_st.k_2);
        config.ff.l_st.X_2 = single(config.ff.l_st.X_2);
        config.ff.l_st.Y_2 = single(config.ff.l_st.Y_2);
        config.ff.l_st.Z_2 = single(config.ff.l_st.Z_2);
        config.ff.l_st.DIR_X_2 = single(config.ff.l_st.DIR_X_2);
        config.ff.l_st.DIR_Y_2 = single(config.ff.l_st.DIR_Y_2);
        config.ff.l_st.DIR_Z_2 = single(config.ff.l_st.DIR_Z_2);
        config.ff.v_st.X_2 = single(config.ff.v_st.X_2);
        config.ff.v_st.Y_2 = single(config.ff.v_st.Y_2);
        config.ff.v_st.Z_2 = single(config.ff.v_st.Z_2);
        config.ff.v_st.DIR_X_2 = single(config.ff.v_st.DIR_X_2);
        config.ff.v_st.DIR_Y_2 = single(config.ff.v_st.DIR_Y_2);
        config.ff.v_st.DIR_Z_2 = single(config.ff.v_st.DIR_Z_2);
    end
    
    if(config.scatter.type == 2)
        config.scatter.f = single(config.scatter.f);
    end
    
    if(config.scatter.type == 3)
        config.scatter.g = single(config.scatter.g);
    end
    
    if(config.sample.direction_type == 2)
        if(isfield(config.sample,'f0'))
            config.sample.f0 = single(config.sample.f0);
        end
    end
    
    if(config.sample.direction_type == 3)
        if(isfield(config.sample,'g0'))
            config.sample.g0 = single(config.sample.g0);
        end
    end
end

%% Scattering function
config.scatter.type = uint32(config.scatter.type);
if(config.scatter.type == 2)
    config.scatter.f = complex(config.scatter.f);
end

end