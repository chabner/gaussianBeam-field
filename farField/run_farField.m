function [u,us] = run_farField(config)

    if(isfield(config,'rng'))
        rng(config.rng);
    end

    if(config.noAttDir)
        if(isfield(config,'nf'))
            [u,us] = MCfieldOnWave( ...
                1/config.MFP,                         ... sigt
                1,                                    ... albedo
                config.box_min,                       ... box_min
                config.box_max,                       ... box_max
                config.iterationsRender,              ... maxItr
                config.multiplePaths,                 ... pathsNum
                config.evalWavelenght,                ... lambda
                config.sampleFlag,                    ... smpFlg
                config.sct_type,                      ... sct_type
                config.ampfunc,                       ... ampfunc
                config.ampfunc0,                      ... ampfunc0
                config.l,                             ... l
                config.v,                             ... v
                config.l,                             ... dirl
                config.v,                             ... dirv
                config.nf                             ... refocus
                );
        else
            [u,us] = MCfieldOnWave( ...
                1/config.MFP,                         ... sigt
                1,                                    ... albedo
                config.box_min,                       ... box_min
                config.box_max,                       ... box_max
                config.iterationsRender,              ... maxItr
                config.multiplePaths,                 ... pathsNum
                config.evalWavelenght,                ... lambda
                config.sampleFlag,                    ... smpFlg
                config.sct_type,                      ... sct_type
                config.ampfunc,                       ... ampfunc
                config.ampfunc0,                      ... ampfunc0
                config.l,                             ... l
                config.v                              ... v
                );
        end
    else
        if(isfield(config,'nf'))
            [u,us] = MCfieldOnWave( ...
                1/config.MFP,                         ... sigt
                1,                                    ... albedo
                config.box_min,                       ... box_min
                config.box_max,                       ... box_max
                config.iterationsRender,              ... maxItr
                config.multiplePaths,                 ... pathsNum
                config.evalWavelenght,                ... lambda
                config.sampleFlag,                    ... smpFlg
                config.sct_type,                      ... sct_type
                config.ampfunc,                       ... ampfunc
                config.ampfunc0,                      ... ampfunc0
                config.l,                             ... l
                config.v,                             ... v
                config.dir_l,                         ... dirl
                config.dir_v,                         ... dirv
                config.nf                             ... refocus
                );
        else
            [u,us] = MCfieldOnWave( ...
                1/config.MFP,                         ... sigt
                1,                                    ... albedo
                config.box_min,                       ... box_min
                config.box_max,                       ... box_max
                config.iterationsRender,              ... maxItr
                config.multiplePaths,                 ... pathsNum
                config.evalWavelenght,                ... lambda
                config.sampleFlag,                    ... smpFlg
                config.sct_type,                      ... sct_type
                config.ampfunc,                       ... ampfunc
                config.ampfunc0,                      ... ampfunc0
                config.l,                             ... l
                config.v,                             ... v
                config.dir_l,                         ... dirl
                config.dir_v                          ... dirv
                );
        end
    end
    

end