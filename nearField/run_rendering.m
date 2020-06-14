function [u,us,xRep,x,w,pxpwVec] = run_rendering(config) 

if(isfield(config,'rng'))
	rng(config.rng);
end

[u,us,~,xRep,~,x,w,pxpwVec] = MCfieldMOvMF( ...
  1/config.MFP,                          ... sigt
  1,                                     ... albedo
  config.box_min,                        ... box_min
  config.box_max,                        ... box_max
  config.iterationsRender,               ... maxItr
  config.multiplePaths,                  ... pathsNum
  config.wavelenght,                     ... lambda
  config.sampleFlag,                     ... smpFlg
  config.smpPreprocess,                  ... smpFunc
  config.movmf,                          ... movmf
  config.sctType,                        ... sct_type
  config.ampfunc,                        ... ampfunc
  config.useGpu,                         ... gpuEnable
  config.apertureVmf_l,                  ... apertureVmf_l
  config.apertureVmf_v                   ... apertureVmf_v
);

end