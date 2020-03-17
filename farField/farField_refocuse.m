function [u_nf,us_nf] = farField_refocuse(config,l_base,v_base) 
%% Build ff directions grid

Nl = size(config.focalPointsL.vector,2);
Nv = numel(config.focalPointsV.base);

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

%% preprocess

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

%% GPU config
if(config.mcGpuL || config.mcGpuV)
    gpuFunc.active = true;
    gridSize = N_v;
    
    gpuFunc.ff_ms = parallel.gpu.CUDAKernel('ff_ms.ptx','ff_ms.cu');
    gpuFunc.ff_ms.GridSize = [ceil(gridSize/gpuFunc.ff_ms.MaxThreadsPerBlock) 1 1];
    gpuFunc.ff_ms.ThreadBlockSize = [gpuFunc.ff_ms.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.ff_ms,'vSize',int32(gridSize));
    setConstantMemory(gpuFunc.ff_ms,'lSize',int32(N_l));
    setConstantMemory(gpuFunc.ff_ms,'box_min',config.box_min);
    setConstantMemory(gpuFunc.ff_ms,'box_max',config.box_max);
    setConstantMemory(gpuFunc.ff_ms,'fw',config.forwardWeight);
    setConstantMemory(gpuFunc.ff_ms,'sigt',1/(2*config.attMFP));
    setConstantMemory(gpuFunc.ff_ms,'g_up',(1 - config.g * config.g)/(4*pi));
    setConstantMemory(gpuFunc.ff_ms,'g_down1',1 + config.g * config.g);
    setConstantMemory(gpuFunc.ff_ms,'g_down2',-2 * config.g);
    
    gpuFunc.ff_single = parallel.gpu.CUDAKernel('ff_single.ptx','ff_single.cu');
    gpuFunc.ff_single.GridSize = [ceil(gridSize/gpuFunc.ff_single.MaxThreadsPerBlock) 1 1];
    gpuFunc.ff_single.ThreadBlockSize = [gpuFunc.ff_single.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.ff_single,'vSize',int32(gridSize));
    setConstantMemory(gpuFunc.ff_single,'lSize',int32(N_l));
    setConstantMemory(gpuFunc.ff_single,'wSize',int32(1));
    setConstantMemory(gpuFunc.ff_single,'box_min',config.box_min);
    setConstantMemory(gpuFunc.ff_single,'box_max',config.box_max);
    setConstantMemory(gpuFunc.ff_single,'sigt',1/(2*config.attMFP));
else
    gpuFunc.active = false;
end

%% refocus config

focalPointsL = config.focalPointsL.vector;
focalPointsV = config.focalPointsV.vector;

kappaL = 1/config.mask_varL^2;
kappaV = 1/config.mask_varV^2;

focalDirectionL = config.focalDirectionsL.vector;
focalDirectionV = config.focalDirectionsV.vector;

w_l = (kappaL/(2*pi)) .* exp( -kappaL + ...
    kappaL * (l(1,:) .* focalDirectionL(1,:)' + l(2,:) .* focalDirectionL(2,:)' + l(3,:) .* focalDirectionL(3,:)') + ...
    1i * 2*pi * (l(1,:) .* focalPointsL(1,:)' + l(2,:) .* focalPointsL(2,:)' + l(3,:) .* focalPointsL(3,:)')) .* ...
    sqrt(1 - (l(1,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1) - (l(2,:).^2) ./ (l(1,:).^2 + l(2,:).^2 - 1));

if(gpuFunc.active)
    w_l = gpuArray(complex(w_l));
    af_ang_vl = gpuArray(af_ang_vl);
    
    gpuFunc.v.x = gpuArray(v(1,:));
    gpuFunc.v.y = gpuArray(v(2,:));
    gpuFunc.v.z = gpuArray(v(3,:));
end

%% run code
for lightNum = 1:1:Nl
    
    if(isfield(config,'rng'))
        rng(config.rng);
    end
    
    if(gpuFunc.active)
        gpuFunc.dirv.x = repmat(gpuArray(config.focalDirectionsV.vector(1,lightNum)),size(gpuFunc.v.x));
        gpuFunc.dirv.y = repmat(gpuArray(config.focalDirectionsV.vector(2,lightNum)),size(gpuFunc.v.y));
        gpuFunc.dirv.z = repmat(gpuArray(config.focalDirectionsV.vector(3,lightNum)),size(gpuFunc.v.z));
    end
    
    [u_ff,us_ff] = MCfieldOnWave( ...
        gpuFunc, ...
        af_ang_vl, ...
        conj(w_l(lightNum,:)) , ...
        config.focalDirectionsL.vector(:,lightNum),    ... dirl
        config.focalDirectionsV.vector(:,lightNum),    ... dirv
        [1/config.scattgMFP,1/config.attMFP] ,  ... sigt
        1,                                      ... albedo
        config.box_min,                         ... box_min
        config.box_max,                         ... box_max
        l,                                      ... l
        v,                                      ... v
        config.iterationsRender,                ... maxItr
        config.wavelenght,                      ... lambda
        config.sampleFlag,                      ... smpFlg
        config.sct_type,                        ... sct_type
        config.ampfunc,                         ... ampfunc
        config.ampfunc0                         ... ampfunc0
    );
    
    w_v = (kappaV/(2*pi)) .* exp( -kappaV + ...
        kappaV * (v(1,:) .* focalDirectionV(1,lightNum) + v(2,:) .* focalDirectionV(2,lightNum) + v(3,:) .* focalDirectionV(3,lightNum)) + ...
        1i * 2*pi * (v(1,:) .* focalPointsV(1,:)' + v(2,:) .* focalPointsV(2,:)' + v(3,:) .* focalPointsV(3,:)')) .* ...
        sqrt(1 - (v(1,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1) - (v(2,:).^2) ./ (v(1,:).^2 + v(2,:).^2 - 1));
    
    u_nf(:,:,lightNum) = reshape(gather(w_v * u_ff * dl * dv), Nv, Nv);
    us_nf(:,:,lightNum) = reshape(gather(w_v * us_ff * dl * dv), Nv, Nv);
end


end