function [u,us,xRep,xVec,wVec,pxpwVec] = run_rendering2D(config) 

if(isfield(config,'rng'))
	rng(config.rng);
end

if(isfield(config,'apertureVmf_l_2'))
    [u,us,~,xRep,xVec,wVec,pxpwVec] = MCfieldMOvMF2D( ...
      [1/config.scattgMFP,1/config.attMFP] , ... sigt
      config.albedo,                         ... albedo
      config.box_min,                        ... box_min
      config.box_max,                        ... box_max
      config.iterationsRender,               ... maxItr
      config.wavelenght,                     ... lambda
      config.movmf,                          ... movmf
      config.sctType,                        ... sct_type
      config.ampfunc,                        ... ampfunc
      config.sampleFlag,                     ... smpFlg
      config.smpPreprocess,                  ... smpFunc
      config.apertureVmf_l,                  ... apertureVmf_l_1
      config.focalDirectionsL.vector,        ... dirl_1
      config.apertureVmf_v,                  ... apertureVmf_v_1
      config.focalDirectionsV.vector,        ... dirv_1
      config.apertureVmf_l_2,                ... apertureVmf_l_2
      config.focalDirectionsL.vector_2,      ... dirl_2
      config.apertureVmf_v_2,                ... apertureVmf_v_2
      config.focalDirectionsV.vector_2       ... dirv_2
    );
else
    [u,us,~,xRep,xVec,wVec,pxpwVec] = MCfieldMOvMF2D( ...
      [1/config.scattgMFP,1/config.attMFP] , ... sigt
      config.albedo,                         ... albedo
      config.box_min,                        ... box_min
      config.box_max,                        ... box_max
      config.iterationsRender,               ... maxItr
      config.wavelenght,                     ... lambda
      config.movmf,                          ... movmf
      config.sctType,                        ... sct_type
      config.ampfunc,                        ... ampfunc
      config.sampleFlag,                     ... smpFlg
      config.smpPreprocess,                  ... smpFunc
      config.apertureVmf_l,                  ... apertureVmf_l_1
      config.focalDirectionsL.vector,        ... dirl_1
      config.apertureVmf_v,                  ... apertureVmf_v_1
      config.focalDirectionsV.vector         ... dirv_1
    );

end
