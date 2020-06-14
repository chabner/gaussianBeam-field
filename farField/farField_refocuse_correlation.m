function [C_nf,Cs_nf] = farField_refocuse_correlation(config,l_base,v_base,focalPointsL_1,focalPointsL_2, ...
    focalPointsV,zeroPer,direction,focalDirectionL,focalDirectionV,maxApertureVal,binaryAperture)
%% Build ff directions grid

if(nargin < 7)
    zeroPer = 0;
end

if(nargin < 8)
    direction = 'full';
end

if(nargin < 9)
    focalDirectionL = [0;0;1];
end

if(nargin < 10)
    focalDirectionV = [0;0;1];
end

if(nargin < 11)
    maxApertureVal = 1;
end

if(nargin < 12)
    binaryAperture = false;
end

Nl = size(focalPointsL_1,2);
Nv = sqrt(size(focalPointsV,2));

dl = l_base(2) - l_base(1);
dl = dl^2;

dv = v_base(2) - v_base(1);
dv = dv^2;

[l_x,l_y] = ndgrid(l_base,l_base);
[v_x,v_y] = ndgrid(v_base,v_base);

l_x = l_x(:)'; l_y = l_y(:)';
v_x = v_x(:)'; v_y = v_y(:)';

badIdx = ((sqrt(l_x.^2 + l_y.^2)) >= maxApertureVal);
l_x = l_x(~badIdx); l_y = l_y(~badIdx);

badIdx = ((sqrt(v_x.^2 + v_y.^2)) >= maxApertureVal);
v_x = v_x(~badIdx); v_y = v_y(~badIdx);

if(strcmp(direction,'full'))
    l = [[l_x;l_y;sqrt(1-l_x.^2-l_y.^2)],[l_x;l_y;-sqrt(1-l_x.^2-l_y.^2)]];
    v = [[v_x;v_y;sqrt(1-v_x.^2-v_y.^2)],[v_x;v_y;-sqrt(1-v_x.^2-v_y.^2)]];
end

if(strcmp(direction,'forward'))
    l = [l_x;l_y;sqrt(1-l_x.^2-l_y.^2)];
    v = [v_x;v_y;sqrt(1-v_x.^2-v_y.^2)];
end

if(strcmp(direction,'backward'))
    l = [l_x;l_y;sqrt(1-l_x.^2-l_y.^2)];
    v = [v_x;v_y;-sqrt(1-v_x.^2-v_y.^2)];
end

%% refocus config

kappaL = 1/config.mask_varL^2;
kappaV = 1/config.mask_varV^2;

w_l_1 = zeros(Nl, size(l,2));
w_l_2 = zeros(Nl, size(l,2));

w_v_1 = zeros(Nl * Nv^2, size(v,2));
w_v_2 = zeros(Nl * Nv^2, size(v,2));

%% run code
for lightNum = 1:1:Nl
    
    L1 = focalPointsL_1(:,lightNum);
    L2 = focalPointsL_2(:,lightNum);
    
    V1 = focalPointsV;
    V2 = focalPointsV;
    
    V1(1:2,:,:) = V1(1:2,:,:) + L1(1:2);
    V2(1:2,:,:) = V2(1:2,:,:) + L2(1:2);
        
    if(~binaryAperture)
        w_l_1(lightNum,:) = (kappaL/(2*pi)) .* exp( -kappaL + ...
            kappaL * (l(1,:) *focalDirectionL(1) + l(2,:) * focalDirectionL(2) + l(3,:) * focalDirectionL(3)) + ...
            1i * 2*pi * (l(1,:) .* L1(1) + l(2,:) .* L1(2) + l(3,:) .* L1(3))) .* ...
            sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));

        w_l_2(lightNum,:) = (kappaL/(2*pi)) .* exp( -kappaL + ...
            kappaL * (l(1,:) *focalDirectionL(1) + l(2,:) * focalDirectionL(2) + l(3,:) * focalDirectionL(3)) + ...
            1i * 2*pi * (l(1,:) .* L2(1) + l(2,:) .* L2(2) + l(3,:) .* L2(3))) .* ...
            sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));

        w_v_1((1 + Nv^2 * (lightNum - 1)):(Nv^2 * lightNum),:) = (kappaV/(2*pi)) .* exp( -kappaV + ...
            kappaV * (v(1,:) *focalDirectionV(1) + v(2,:) * focalDirectionV(2) + v(3,:) * focalDirectionV(3)) + ...
            1i * 2*pi * (v(1,:) .* V1(1,:)' + v(2,:) .* V1(2,:)' + v(3,:) .* V1(3,:)')) .* ...
            sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));

        w_v_2((1 + Nv^2 * (lightNum - 1)):(Nv^2 * lightNum),:) = (kappaV/(2*pi)) .* exp( -kappaV + ...
            kappaV * (v(1,:) *focalDirectionV(1) + v(2,:) * focalDirectionV(2) + v(3,:) * focalDirectionV(3)) + ...
            1i * 2*pi * (v(1,:) .* V2(1,:)' + v(2,:) .* V2(2,:)' + v(3,:) .* V2(3,:)')) .* ...
            sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));
    else
        
        w_l_1(lightNum,:) = exp( ...
            1i * 2*pi * (l(1,:) .* L1(1) + l(2,:) .* L1(2) + l(3,:) .* L1(3))) .* ...
            sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));

        w_l_2(lightNum,:) = exp(...
            1i * 2*pi * (l(1,:) .* L2(1) + l(2,:) .* L2(2) + l(3,:) .* L2(3))) .* ...
            sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));

        w_v_1((1 + Nv^2 * (lightNum - 1)):(Nv^2 * lightNum),:) = exp(...
            1i * 2*pi * (v(1,:) .* V1(1,:)' + v(2,:) .* V1(2,:)' + v(3,:) .* V1(3,:)')) .* ...
            sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));

        w_v_2((1 + Nv^2 * (lightNum - 1)):(Nv^2 * lightNum),:) = exp( ...
            1i * 2*pi * (v(1,:) .* V2(1,:)' + v(2,:) .* V2(2,:)' + v(3,:) .* V2(3,:)')) .* ...
            sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));
    end
end

if(isfield(config,'rng'))
    rng(config.rng);
end

dirl = repmat(focalDirectionL,[1,size(l,2)]);
dirv = repmat(focalDirectionV,[1,size(v,2)]);
    
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

C_nf = zeros(Nv,Nv,Nl);
Cs_nf = zeros(Nv,Nv,Nl);

for lightNum = 1:1:Nl
    curr_u_ff = u_ff((1 + Nv^2 * (lightNum - 1)):(Nv^2 * lightNum),lightNum);
    curr_us_ff = us_ff((1 + Nv^2 * (lightNum - 1)):(Nv^2 * lightNum),lightNum);
    C_nf(:,:,lightNum) = reshape(curr_u_ff * (dl^2) * (dv^2), Nv, Nv);
    Cs_nf(:,:,lightNum) = reshape(curr_us_ff * (dl^2) * (dv^2), Nv, Nv);
end

end

