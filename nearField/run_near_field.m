function [u,us] = run_near_field(config) 

if(numel(config.simulation.gpuNum) == 1)
    [u,us,normalizationFactor,totalIter] = near_field_function( ...
        config.simulation,       ... simulation struct
        config.medium,           ... medium struct
        config.aperture,         ... aperture struct
        config.sample,           ... sampling struct
        config.scatter,          ... scattering struct
        config.movmf,            ... mixture of vmf of the scattering function
        config.nf.wavenumber_st, ... k number
        config.nf.l_st,          ... illumination
        config.nf.v_st           ... view
    );
else  
    gpuList = config.simulation.gpuNum;
    curr_u = cell(1,numel(gpuList));
    curr_us = cell(1,numel(gpuList));
    curr_totalIter = cell(1,numel(gpuList));
    curr_normalizationFactor = cell(1,numel(gpuList));
    
    parfor iter_num = 1:numel(config.simulation.gpuNum)
        pfor_config = config;
        pfor_config.simulation.gpuNum = gpuList(iter_num);
        
        [curr_u{iter_num},curr_us{iter_num},curr_normalizationFactor{iter_num},curr_totalIter{iter_num}] = near_field_function( ...
            pfor_config.simulation,       ... simulation struct
            pfor_config.medium,           ... medium struct
            pfor_config.aperture,         ... aperture struct
            pfor_config.sample,           ... sampling struct
            pfor_config.scatter,          ... scattering struct
            pfor_config.movmf,            ... mixture of vmf of the scattering function
            pfor_config.nf.wavenumber_st, ... k number
            pfor_config.nf.l_st,          ... illumination
            pfor_config.nf.v_st           ... view
    );
    end
    
    u = 0;
    us = 0;
    totalIter = 0;
    normalizationFactor = curr_normalizationFactor{1};
    
    for iter_num = 1:numel(config.gpuNum)
        u = u + curr_u{iter_num};
        us = us + curr_us{iter_num};
        totalIter = totalIter + curr_totalIter{iter_num};
    end
end

correlationFlag = isfield(config.nf.l_st,'P1_2');

if(correlationFlag)
    u = u * normalizationFactor/totalIter;
    us = us * (normalizationFactor/totalIter);
else
    u = u * sqrt(normalizationFactor/totalIter);
    us = us * sqrt(normalizationFactor/totalIter);
end

end