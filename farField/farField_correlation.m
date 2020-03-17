function [C_ff,Cs_ff] = farField_correlation(config,l_base,v_base,itersNum) 
%% GPU code
if(config.mcGpuL || config.mcGpuV)
    gpuFunc.active = true;
    gridSize = 2*(numel(v_base)^2);
    
    gpuFunc.ff_ms = parallel.gpu.CUDAKernel('ff_ms.ptx','ff_ms.cu');
    gpuFunc.ff_ms.GridSize = [ceil((2*gridSize)/gpuFunc.ff_ms.MaxThreadsPerBlock) 1 1];
    gpuFunc.ff_ms.ThreadBlockSize = [gpuFunc.ff_ms.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.ff_ms,'vSize',int32(gridSize));
    setConstantMemory(gpuFunc.ff_ms,'lSize',int32(2));
    setConstantMemory(gpuFunc.ff_ms,'box_min',config.box_min);
    setConstantMemory(gpuFunc.ff_ms,'box_max',config.box_max);
    setConstantMemory(gpuFunc.ff_ms,'fw',config.forwardWeight);
    setConstantMemory(gpuFunc.ff_ms,'sigt',1/(2*config.attMFP));
    setConstantMemory(gpuFunc.ff_ms,'g_up',(1 - config.g * config.g)/(4*pi));
    setConstantMemory(gpuFunc.ff_ms,'g_down1',1 + config.g * config.g);
    setConstantMemory(gpuFunc.ff_ms,'g_down2',-2 * config.g);
    
    gpuFunc.ff_single = parallel.gpu.CUDAKernel('ff_single.ptx','ff_single.cu');
    gpuFunc.ff_single.GridSize = [ceil((2*gridSize)/gpuFunc.ff_single.MaxThreadsPerBlock) 1 1];
    gpuFunc.ff_single.ThreadBlockSize = [gpuFunc.ff_single.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.ff_single,'vSize',int32(gridSize));
    setConstantMemory(gpuFunc.ff_single,'lSize',int32(2));
    setConstantMemory(gpuFunc.ff_single,'wSize',int32(2));
    setConstantMemory(gpuFunc.ff_single,'box_min',config.box_min);
    setConstantMemory(gpuFunc.ff_single,'box_max',config.box_max);
    setConstantMemory(gpuFunc.ff_single,'sigt',1/(2*config.attMFP));
else
    gpuFunc.active = false;
end

%% GB code
[v_x,v_y] = ndgrid(v_base,v_base);
v_x = v_x(:).'; v_y = v_y(:).';

N_l = numel(l_base);
N_v = numel(v_base).^2;

C_ff = zeros(N_v,N_l);
Cs_ff = C_ff;

if(gpuFunc.active)
    C_ff = gpuArray(complex(C_ff));
    Cs_ff = gpuArray(complex(Cs_ff));
end

for lNum = 1:1:numel(l_base)
    l_x = l_base(lNum)/2;
    l = [-l_x,l_x;0,0;sqrt(1-l_x.^2),sqrt(1-l_x.^2)];
%     l = [0,l_x;0,0;1,sqrt(1-l_x.^2)];
    
    v1_x = sin((asin(v_x) - asin(l_x)));
    v2_x = sin((asin(v_x) + asin(l_x)));
    
    v1 = [v1_x;v_y;sqrt(1-v1_x.^2-v_y.^2)];
    v2 = [v2_x;v_y;sqrt(1-v2_x.^2-v_y.^2)];
    
%     v1 = [(v_x);v_y;sqrt(1-(v_x).^2-v_y.^2)];
%     v2 = v1 + l(:,2) - l(:,1);
%     v2(3,:) = sqrt(1 - v2(1,:).^2 - v2(2,:).^2);
    v = [v1,v2];
    
    
    
    af_ang_vl = zeros(2*N_v, 2);
    if config.sct_type>1
        for j=1:2
            af_ang_vl(:,j)=evalampfunc_general(l(:,j)'*v,config.sct_type,config.ampfunc,config.dimNum);
        end
    else
        af_ang_vl=evalampfunc_general(0,config.sct_type,config.ampfunc,config.dimNum);
    end
    
    Wl = eye(2);
    
    if(gpuFunc.active)
        gpuFunc.v.x = gpuArray(v(1,:));
        gpuFunc.v.y = gpuArray(v(2,:));
        gpuFunc.v.z = gpuArray(v(3,:));

        gpuFunc.dirv.x = gpuArray(v(1,:));
        gpuFunc.dirv.y = gpuArray(v(2,:));
        gpuFunc.dirv.z = gpuArray(v(3,:));
        
        Wl = gpuArray(complex(Wl));
        af_ang_vl = gpuArray(af_ang_vl);
    end

    for iterNum = 1:1:itersNum
        [u_ff,us_ff] = MCfieldOnWave( ...
            gpuFunc , ...
            af_ang_vl, ...
            Wl , ...
            l,    ... dirl
            v,    ... dirv
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
    
        u_ff1 = u_ff(1:N_v,1);
        u_ff2 = u_ff(N_v+1:end,2);
        C_ff(:,lNum) = C_ff(:,lNum) + u_ff1 .* conj(u_ff2);
        
        us_ff1 = us_ff(1:N_v,1);
        us_ff2 = us_ff(N_v+1:end,2);
        Cs_ff(:,lNum) = Cs_ff(:,lNum) + us_ff1 .* conj(us_ff2); 
    end
end

C_ff = gather(C_ff);
Cs_ff = gather(Cs_ff);

C_ff = reshape(C_ff,[sqrt(N_v),sqrt(N_v),N_l]);
Cs_ff = reshape(Cs_ff,[sqrt(N_v),sqrt(N_v),N_l]);

end