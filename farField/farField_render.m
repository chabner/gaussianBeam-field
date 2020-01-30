function [u_nf,us_nf] = farField_render(config,l_base,v_base) 
%% GB code
% config = buildConfig(config);

Nl = size(config.focalPointsL,2);
Nv = numel(config.focalPointsV_base);

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

badIdx = (l_x.^2 + l_y.^2 > 0.7);
l_x = l_x(~badIdx); l_y = l_y(~badIdx);

badIdx = (v_x.^2 + v_y.^2 > 0.7);
v_x = v_x(~badIdx); v_y = v_y(~badIdx);

l = [l_x;l_y;sqrt(1-l_x.^2-l_y.^2)];
v = [v_x;v_y;sqrt(1-v_x.^2-v_y.^2)];

N_l = size(l,2);
N_v = size(v,2);
af_ang_vl = zeros(N_v, N_l);
if config.sct_type>1
    for j=1:N_v
        af_ang_vl(:,j)=evalampfunc_general(l(:,j)'*v,config.sct_type,config.ampfunc,config.dimNum);
    end
else
    af_ang_vl=evalampfunc_general(0,config.sct_type,config.ampfunc,config.dimNum);
end

focalPointsL = config.focalPointsL;
focalPointsV = config.focalPointsV;

kappaL = 1/config.mask_varL^2;
kappaV = 1/config.mask_varV^2;

focalDirectionL = config.focalDirectionsL;
focalDirectionV = config.focalDirectionsV;

for lightNum = 1:1:Nl
    
    if(isfield(config,'rng'))
        rng(config.rng);
    end
    
    w_l = (kappaL/(2*pi)) .* exp( -kappaL + ...
        kappaL * (l(1,:) .* focalDirectionL(1,lightNum) + l(2,:) .* focalDirectionL(2,lightNum) + l(3,:) .* focalDirectionL(3,lightNum)) + ...
        1i * 2*pi * (l(1,:) .* focalPointsL(1,lightNum)' + l(2,:) .* focalPointsL(2,lightNum)' + l(3,:) .* focalPointsL(3,lightNum)')) .* ...
        sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));
    
    w_v = (kappaV/(2*pi)) .* exp( -kappaV + ...
        kappaV * (v(1,:) .* focalDirectionV(1,lightNum) + v(2,:) .* focalDirectionV(2,lightNum) + v(3,:) .* focalDirectionV(3,lightNum)) + ...
        1i * 2*pi * (v(1,:) .* focalPointsV(1,:)' + v(2,:) .* focalPointsV(2,:)' + v(3,:) .* focalPointsV(3,:)')) .* ...
        sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));
    
    [u_ff,us_ff] = MCfieldOnWave( ...
        af_ang_vl, ...
        conj(w_l) , ...
        config.focalDirectionsL(:,lightNum),    ... dirl
        config.focalDirectionsV(:,lightNum),    ... dirv
        [1/config.scattgMFP,1/config.attMFP] ,  ... sigt
        1,                                      ... albedo
        config.box_min,                         ... box_min
        config.box_max,                         ... box_max
        l,                                      ... l
        v,                                      ... v
        true,                                   ... is_ff_l
        true,                                   ... is_ff_v
        config.iterationsRender,                ... maxItr
        config.wavelenght,                      ... lambda
        false,                                  ... doCBS
        config.sampleFlag,                      ... smpFlg
        config.sct_type,                        ... sct_type
        config.ampfunc                          ... ampfunc
    );
    
    u_nf(:,:,lightNum) = reshape(w_v * u_ff * dl * dv, Nv, Nv);
    us_nf(:,:,lightNum) = reshape(w_v * us_ff * dl * dv, Nv, Nv);
end


end