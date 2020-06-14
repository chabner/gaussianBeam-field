function [u_nf,us_nf] = farField_refocuse(config,l_base,v_base,focalPointsL,focalPointsV,focalDirectionL,focalDirectionV,zParam) 
if(nargin == 7)
    zParam = 'full';
end
%% Build ff directions grid

Nl = size(focalPointsL,2);
Nv = sqrt(size(focalPointsV,2));

u_nf = zeros(Nv,Nv,Nl);
us_nf = zeros(Nv,Nv,Nl);

dl = l_base(2) - l_base(1);
dl = dl^2;

dv = v_base(2) - v_base(1);
dv = dv^2;

[l_x,l_y] = ndgrid(l_base,l_base);
[v_x,v_y] = ndgrid(v_base,v_base);

l_x = l_x(:)'; l_y = l_y(:)';
v_x = v_x(:)'; v_y = v_y(:)';

badIdx = (l_x.^2 + l_y.^2 >= 1);
l_x = l_x(~badIdx); l_y = l_y(~badIdx);

badIdx = (v_x.^2 + v_y.^2 >= 1);
v_x = v_x(~badIdx); v_y = v_y(~badIdx);

l = [[l_x;l_y;sqrt(1-l_x.^2-l_y.^2)],[l_x;l_y;-sqrt(1-l_x.^2-l_y.^2)]];

if(strcmp(zParam,'full'))
   v = [[v_x;v_y;sqrt(1-v_x.^2-v_y.^2)],[l_x;l_y;-sqrt(1-l_x.^2-l_y.^2)]]; 
end

if(strcmp(zParam,'forward'))
   v = [v_x;v_y;sqrt(1-v_x.^2-v_y.^2)]; 
end

if(strcmp(zParam,'backward'))
   v = [v_x;v_y;-sqrt(1-v_x.^2-v_y.^2)]; 
end

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
    
    [u_ff,us_ff] = MCfieldOnWave( ...
        1,                                    ... gpuNum
        1/config.MFP,                         ... sigt
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
        []                                    ... wl_1
        );
    
    w_l = (kappaL/(2*pi)) .* exp( -kappaL + ...
        kappaL * (l(1,:) .* focalDirectionL(1,lightNum)' + l(2,:) .* focalDirectionL(2,lightNum)' + l(3,:) .* focalDirectionL(3,lightNum)') + ...
        1i * 2*pi * (l(1,:) .* focalPointsL(1,lightNum)' + l(2,:) .* focalPointsL(2,lightNum)' + l(3,:) .* focalPointsL(3,lightNum)')) .* ...
        sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));
    
    w_v = (kappaV/(2*pi)) .* exp( -kappaV + ...
        kappaV * (v(1,:) .* focalDirectionV(1,lightNum) + v(2,:) .* focalDirectionV(2,lightNum) + v(3,:) .* focalDirectionV(3,lightNum)) + ...
        1i * 2*pi * (v(1,:) .* focalPointsV(1,:)' + v(2,:) .* focalPointsV(2,:)' + v(3,:) .* focalPointsV(3,:)')) .* ...
        sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));
    
    u_nf(:,:,lightNum) = reshape(w_v * u_ff * w_l' * dl * dv, Nv, Nv);
    us_nf(:,:,lightNum) = reshape(w_v * us_ff * w_l' * dl * dv, Nv, Nv);
end

end