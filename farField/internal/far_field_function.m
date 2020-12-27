function [u,us,normalizationFactor,totalIter] = far_field_function(simulation,medium,sampling,scattering,wavenumber,illumination,view)

[u,us,normalizationFactor,totalIter] = ff_mex( ...
    simulation,     ... simulation struct
    medium,         ... medium struct
    sampling,       ... sampling struct
    scattering,     ... scattering struct
    wavenumber,     ... k number
    illumination,   ... illumination
    view            ... view
);

u = double(u);
us = double(us);
u = u + us;

normalizationFactor = double(normalizationFactor);
totalIter = double(totalIter);

end
