function [u,us,normalizationFactor,totalIter] = near_field_function(simulation,medium,aperture,sampling,scattering,movmf,wavenumber,illumination,view)

[u,us,normalizationFactor,totalIter,min_px0,min_pw0] = nf_mex( ...
    simulation,     ... simulation struct
    medium,         ... medium struct
    aperture,       ... aperture struct
    sampling,       ... sampling struct
    scattering,     ... scattering struct
    movmf,          ... mixture of vmf of the scattering function
    wavenumber,     ... k number
    illumination,   ... illumination
    view            ... view
);

% min_px0
% min_pw0

u = double(u);
us = double(us);
u = u + us;

normalizationFactor = double(normalizationFactor);
totalIter = double(totalIter);


end
