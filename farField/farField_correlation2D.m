function [C_ff,Cs_ff] = farField_correlation2D(config,theta_l,theta_v) 
%% Build ff directions grid

l = [sin(theta_l);cos(theta_l)];
v = [sin(theta_v);cos(theta_v)];

%% run code
if(isfield(config,'rng'))
    rng(config.rng);
end
    
[C_ff,Cs_ff] = MCfieldOnWave( ...
    1,                                    ... gpuNum
    [1/config.scattgMFP,1/config.attMFP], ... sigt
    1,                                    ... albedo
    config.box_min,                       ... box_min
    config.box_max,                       ... box_max
    config.iterationsRender,              ... maxItr
    config.wavelenght,                    ... lambda
    config.sampleFlag,                    ... smpFlg
    config.sct_type,                      ... sct_type
    config.ampfunc,                       ... ampfunc
    config.ampfunc0,                      ... ampfunc0
    l,                                    ... l_1
    v,                                    ... v_1
    l,                                 ... dirl_1
    v,                                 ... dirv_1
    [],                             ... wl_1
    [], ...
    l,                                    ... l_1
    v,                                    ... v_1
    l,                                 ... dirl_1
    v,                                 ... dirv_1
    [],                             ... wl_1
    [] ...
    );

end

