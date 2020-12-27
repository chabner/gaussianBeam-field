function [u,us] = run_far_field(config) 

if(numel(config.simulation.gpuNum) == 1)
    rng('shuffle');
    
    [u,us,normalizationFactor,totalIter] = far_field_function( ...
      config.simulation,       ... simulation
      config.medium,           ... medium
      config.sample,           ... sampling
      config.scatter,          ... scattering
      config.ff.wavenumber_st, ... wavenumber
      config.ff.l_st,          ... illumination
      config.ff.v_st           ... view
    );
end

correlationFlag = isfield(config.ff.l_st,'X_2');

if(correlationFlag)
    u = u * normalizationFactor/totalIter;
    us = us * (normalizationFactor/totalIter);
else
    u = u * sqrt(normalizationFactor/totalIter);
    us = us * sqrt(normalizationFactor/totalIter);
end

end