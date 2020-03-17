function [u,us,xRep] = run_rendering_mog(config,gpuFunc) 

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
  config.focalDirectionsL.vector,        ... dirl
  config.focalDirectionsV.vector,        ... dirv
  config.iterationsRender,               ... maxItr
  config.wavelenght,                     ... lambda
  config.sampleFlag,                     ... smpFlg
  config.smpPreprocess,                  ... smpFunc
  config.movmf,                          ... movmf
  config.sct_type,                       ... sct_type
  config.ampfunc,                        ... ampfunc
  gpuFunc                                ... gpuFunc
);

% if(config.mcGpuV || config.mcGpuL)
%     u = gather(u);
%     us = gather(us);
% end

end