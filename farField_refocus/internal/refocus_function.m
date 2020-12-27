function [u,us,normalizationFactor,totalIter] = refocus_function(simulation,medium,aperture,refocus,scattering,wavenumber,illumination,view)

[u,us,normalizationFactor,totalIter] = refocus_mex( ...
    simulation,     ... simulation struct
    medium,         ... medium struct
    aperture,       ... aperture struct
    refocus,        ... refocus struct
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
