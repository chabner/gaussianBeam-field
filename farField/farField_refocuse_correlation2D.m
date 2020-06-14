function [C_nf,Cs_nf] = farField_refocuse_correlation2D(config,theta_l,theta_v,focalPointsL_1,focalPointsL_2,focalPointsV,zeroPer,dirFlag) 
%% Build ff directions grid

if(nargin < 6)
    zeroPer = 0;
end

if(nargin < 7)
    dirFlag = true;
end

Nl = size(focalPointsL_1,2);
Nv = size(focalPointsV,2);

l = [sin(theta_l);cos(theta_l)];
v = [sin(theta_v);cos(theta_v)];

dthetav = theta_v(2) - theta_v(1);
dthetal = theta_l(2) - theta_l(1);

% l = [[l_base;sqrt(1-l_base.^2)],[l_base;-sqrt(1-l_base.^2)]];
% v = [[v_base;sqrt(1-v_base.^2)],[v_base;-sqrt(1-v_base.^2)]];

%% refocus config

kappaL = 1/config.mask_varL^2;
kappaV = 1/config.mask_varV^2;

w_l_1 = zeros(Nl, size(l,2));
w_l_2 = zeros(Nl, size(l,2));

w_v_1 = zeros(Nl * Nv, size(v,2));
w_v_2 = zeros(Nl * Nv, size(v,2));


%% run code
for lightNum = 1:1:Nl
    if(dirFlag)
        dirl = 0 * l;
        dirl(end,:) = 1;
    
        dirv = 0 * v;
        dirv(end,:) = 1;
    else
        dirl = l;
        dirv = v;        
    end


    
    L1 = focalPointsL_1(:,lightNum);
    L2 = focalPointsL_2(:,lightNum);
    
    V1 = focalPointsV;
    V2 = focalPointsV;
    
    V1(1,:) = V1(1,:) + L1(1);
    V2(1,:) = V2(1,:) + L2(1);
    
    
    w_l_1(lightNum,:) = (1/(2*pi*besseli(0,kappaL,1))) .* exp( -kappaL + ...
        kappaL * (l(2,:)) + ...
        1i * 2*pi * (l(1,:) .* L1(1) + l(2,:) .* L1(2)));
    
    w_l_2(lightNum,:) = (1/(2*pi*besseli(0,kappaL,1))) .* exp( -kappaL + ...
        kappaL * (l(2,:)) + ...
        1i * 2*pi * (l(1,:) .* L2(1) + l(2,:) .* L2(2)));
    
    w_v_1((1 + Nv * (lightNum - 1)):(Nv * lightNum),:) = (1/(2*pi*besseli(0,kappaV,1))) .* exp( -kappaV + ...
        kappaV * (v(2,:)) + ...
        1i * 2*pi * (v(1,:) .* V1(1,:)' + v(2,:) .* V1(2,:)'));
    
    w_v_2((1 + Nv * (lightNum - 1)):(Nv * lightNum),:) = (1/(2*pi*besseli(0,kappaV,1))) .* exp( -kappaV + ...
        kappaV * (v(2,:)) + ...
        1i * 2*pi * (v(1,:) .* V2(1,:)' + v(2,:) .* V2(2,:)'));
end
if(isfield(config,'rng'))
    rng(config.rng);
end
    
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
    conj(w_l_1),                             ... wl_1
    w_v_1, ...
    l,                                    ... l_1
    v,                                    ... v_1
    dirl,                                 ... dirl_1
    dirv,                                 ... dirv_1
    conj(w_l_2),                             ... wl_1
    w_v_2, ...
    zeroPer ...
    );
    
C_nf = zeros(Nv,Nl);
Cs_nf = zeros(Nv,Nl);

for lightNum = 1:1:Nl
    curr_u_ff = u_ff((1 + Nv * (lightNum - 1)):(Nv * lightNum),lightNum);
    curr_us_ff = us_ff((1 + Nv * (lightNum - 1)):(Nv * lightNum),lightNum);
    C_nf(:,lightNum) = curr_u_ff * (dthetal * dthetav)^2;
    Cs_nf(:,lightNum) = curr_us_ff * (dthetal * dthetav)^2;
end

end

