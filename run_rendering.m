function [u,us,xRep] = run_rendering(config) 

if(isfield(config,'rng'))
	rng(config.rng);
end

[u,us,~,~,xRep] = MCfieldGaussianBeam( ...
  [1/config.scattgMFP,1/config.attMFP] , ... sigt
  1,                                     ... albedo
  config.box_min,                        ... box_min
  config.box_max,                        ... box_max
  config.apertureVmf_l,                  ... apertureVmf_l
  config.apertureVmf_v,                  ... apertureVmf_v
  1,                                     ... signl
  1,                                     ... signv
  config.mask_varL,                      ... varl
  config.mask_varV,                      ... varv
  config.focalDirectionsL,               ... dirl
  config.focalDirectionsV,               ... dirv
  config.iterationsRender,               ... maxItr
  config.wavelenght,                     ... lambda
  config.sampleFlag,                     ... smpFlg
  config.smpPreprocess,                  ... smpFunc
  config.movmf,                          ... movmf
  config.sct_type,                       ... sct_type
  config.ampfunc,                        ... ampfunc
  config.ampfunc0                        ... ampfunc0
);

if(config.mcGpuV || config.mcGpuL)
    u = gather(u);
    us = gather(us);
end

Nv = numel(config.focalPointsV_base);
Nl = numel(config.focalPointsL_base);

u = reshape(u,Nv,Nv,Nl);
us = reshape(us,Nv,Nv,Nl);


end