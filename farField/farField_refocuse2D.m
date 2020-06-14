function [u_nf,us_nf] = farField_refocuse2D(config,l_base,v_base,focalPointsL,focalPointsV,focalDirectionL,focalDirectionV) 
%% Build ff directions grid

Nl = size(focalPointsL,2);
Nv = size(focalPointsV,2);

u_nf = zeros(Nv,Nl);
us_nf = zeros(Nv,Nl);

dl = l_base(2) - l_base(1);
dv = v_base(2) - v_base(1);

l_base = l_base(:).';
v_base = v_base(:).';

l = [[l_base;sqrt(1-l_base.^2)],[l_base;-sqrt(1-l_base.^2)]];
v = [[v_base;sqrt(1-v_base.^2)],[v_base;-sqrt(1-v_base.^2)]];

%% refocus config

kappaL = 1/config.mask_varL^2;
kappaV = 1/config.mask_varV^2;

%% run code
for lightNum = 1:1:Nl
    
    if(isfield(config,'rng'))
        rng(config.rng);
    end
    
    dirl = repmat(focalDirectionL(:,lightNum) , 1, size(l,2));
    dirv = repmat(focalDirectionV(:,lightNum) , 1, size(v,2));
    
    w_l = (1/(2*pi*besseli(0,kappaL,1))) .* exp( -kappaL + ...
        kappaL * (l(1,:) .* focalDirectionL(1,lightNum)' + l(2,:) .* focalDirectionL(2,lightNum)') + ...
        1i * 2*pi * (l(1,:) .* focalPointsL(1,lightNum)' + l(2,:) .* focalPointsL(2,lightNum)'));
    
    [u_ff,us_ff] = MCfieldOnWave( ...
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
        dirl,                                 ... dirl_1
        dirv,                                 ... dirv_1
        conj(w_l)                             ... wl_1
        );
    
    w_v = (1/(2*pi*besseli(0,kappaV,1))) .* exp( -kappaV + ...
        kappaV * (v(1,:) .* focalDirectionV(1,lightNum) + v(2,:) .* focalDirectionV(2,lightNum)) + ...
        1i * 2*pi * (v(1,:) .* focalPointsV(1,:)' + v(2,:) .* focalPointsV(2,:)'));
    
    u_nf(:,lightNum) = w_v * u_ff * dl * dv;
    us_nf(:,lightNum) = w_v * us_ff * dl * dv;
end

end