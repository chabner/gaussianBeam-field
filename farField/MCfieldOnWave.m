function [u,us]=MCfieldOnWave(gpuNum,sigt,albedo,box_min,box_max,maxItr,lambda,smpFlg,sct_type,ampfunc,ampfunc0, ...
    l_1,v_1,dirl_1,dirv_1,wl_1,wv_1, ...
    l_2,v_2,dirl_2,dirv_2,wl_2,wv_2,zeroPercent)
%%

correlationFlag = (nargin > 17);
noWv = (nargin == 16);
if(noWv)
    wv_1 = 1;
end

if(nargin > 16 && isempty(wv_1))
    wv_1 = 1;
end

if(nargin > 22 && isempty(wv_2))
    wv_2 = 1;
end

noWl = isempty(wl_1);

if(noWl)
    wl_1 = 1;
    wl_2 = 1;
end

if(nargin < 24)
    zeroPercent = 0;
end

dim = size(box_min,1);

%% GPU funcs
reset(gpuDevice(gpuNum));

Nv = size(v_1,2);
Nl = size(l_1,2);
Nw = size(wl_1,1);

if(noWl)
    Nw = Nl;
end

if(dim == 3)
    if(sct_type == 3)
        ff_ms = parallel.gpu.CUDAKernel('ff_ms.ptx','ff_ms.cu');
        ff_ms.GridSize = [ceil((Nv * Nw)/ff_ms.MaxThreadsPerBlock) 1 1];
        ff_ms.ThreadBlockSize = [ff_ms.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_ms,'vSize',int32(Nv));
        setConstantMemory(ff_ms,'lSize',int32(Nw));
        setConstantMemory(ff_ms,'box_min',box_min);
        setConstantMemory(ff_ms,'box_max',box_max);
        setConstantMemory(ff_ms,'fw',ampfunc.forwardWeight);
        setConstantMemory(ff_ms,'sigt',sigt/2);
        setConstantMemory(ff_ms,'g_up',(1 - ampfunc.g * ampfunc.g)/(4*pi));
        setConstantMemory(ff_ms,'g_down1',1 + ampfunc.g * ampfunc.g);
        setConstantMemory(ff_ms,'g_down2',-2 * ampfunc.g);
    end
    
    if(sct_type == 2)
        ff_ms_tabulated = parallel.gpu.CUDAKernel('ff_ms_tabulated.ptx','ff_ms_tabulated.cu');
        ff_ms_tabulated.GridSize = [ceil((Nv * Nw)/ff_ms_tabulated.MaxThreadsPerBlock) 1 1];
        ff_ms_tabulated.ThreadBlockSize = [ff_ms_tabulated.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_ms_tabulated,'vSize',int32(Nv));
        setConstantMemory(ff_ms_tabulated,'lSize',int32(Nw));
        setConstantMemory(ff_ms_tabulated,'nAmpfunc',int32(numel(ampfunc.evalAmp)));
        setConstantMemory(ff_ms_tabulated,'box_min',box_min);
        setConstantMemory(ff_ms_tabulated,'box_max',box_max);
        setConstantMemory(ff_ms_tabulated,'sigt',sigt/2);
    end

    if(~noWl)
        ff_single = parallel.gpu.CUDAKernel('ff_single.ptx','ff_single.cu');
        ff_single.GridSize = [ceil((Nv * Nw)/ff_single.MaxThreadsPerBlock) 1 1];
        ff_single.ThreadBlockSize = [ff_single.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_single,'vSize',int32(Nv));
        setConstantMemory(ff_single,'lSize',int32(Nl));
        setConstantMemory(ff_single,'wSize',int32(Nw));
        setConstantMemory(ff_single,'box_min',box_min);
        setConstantMemory(ff_single,'box_max',box_max);
        setConstantMemory(ff_single,'sigt',sigt/2);
    else
        ff_single_noW = parallel.gpu.CUDAKernel('ff_single_noW.ptx','ff_single_noW.cu');
        ff_single_noW.GridSize = [ceil((Nv * Nw)/ff_single_noW.MaxThreadsPerBlock) 1 1];
        ff_single_noW.ThreadBlockSize = [ff_single_noW.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_single_noW,'vSize',int32(Nv));
        setConstantMemory(ff_single_noW,'lSize',int32(Nl));
        setConstantMemory(ff_single_noW,'box_min',box_min);
        setConstantMemory(ff_single_noW,'box_max',box_max);
        setConstantMemory(ff_single_noW,'sigt',sigt/2);
    end

    if(sct_type == 3)
        af_ang_vl_gpu = parallel.gpu.CUDAKernel('af_ang_vl_gpu.ptx','af_ang_vl_gpu.cu');
        af_ang_vl_gpu.GridSize = [ceil((Nv * Nl)/af_ang_vl_gpu.MaxThreadsPerBlock) 1 1];
        af_ang_vl_gpu.ThreadBlockSize = [af_ang_vl_gpu.MaxThreadsPerBlock 1 1];

        setConstantMemory(af_ang_vl_gpu,'vSize',int32(Nv));
        setConstantMemory(af_ang_vl_gpu,'lSize',int32(Nl));
        setConstantMemory(af_ang_vl_gpu,'fw',ampfunc.forwardWeight);
        setConstantMemory(af_ang_vl_gpu,'g_up',(1 - ampfunc.g * ampfunc.g)/(4*pi));
        setConstantMemory(af_ang_vl_gpu,'g_down1',1 + ampfunc.g * ampfunc.g);
        setConstantMemory(af_ang_vl_gpu,'g_down2',-2 * ampfunc.g);
    end
    
    if(sct_type == 2)
        af_ang_vl_gpu_tabulated = parallel.gpu.CUDAKernel('af_ang_vl_gpu_tabulated.ptx','af_ang_vl_gpu_tabulated.cu');
        af_ang_vl_gpu_tabulated.GridSize = [ceil((Nv * Nl)/af_ang_vl_gpu_tabulated.MaxThreadsPerBlock) 1 1];
        af_ang_vl_gpu_tabulated.ThreadBlockSize = [af_ang_vl_gpu_tabulated.MaxThreadsPerBlock 1 1];

        setConstantMemory(af_ang_vl_gpu_tabulated,'vSize',int32(Nv));
        setConstantMemory(af_ang_vl_gpu_tabulated,'lSize',int32(Nl));
        setConstantMemory(af_ang_vl_gpu_tabulated,'nAmpfunc',int32(numel(ampfunc.evalAmp)));
    end

    wl_1 = gpuArray(complex(wl_1));
    wv_1 = gpuArray(wv_1);
    vx_1 = gpuArray(v_1(1,:));
    vy_1 = gpuArray(v_1(2,:));
    vz_1 = gpuArray(v_1(3,:));
    lx_1 = gpuArray(l_1(1,:));
    ly_1 = gpuArray(l_1(2,:));
    lz_1 = gpuArray(l_1(3,:));
    dirvx_1 = gpuArray(dirv_1(1,:));
    dirvy_1 = gpuArray(dirv_1(2,:));
    dirvz_1 = gpuArray(dirv_1(3,:));

    if(correlationFlag)
        wl_2 = gpuArray(complex(wl_2));
        wv_2 = gpuArray(wv_2);
        vx_2 = gpuArray(v_2(1,:));
        vy_2 = gpuArray(v_2(2,:));
        vz_2 = gpuArray(v_2(3,:));
        lx_2 = gpuArray(l_2(1,:));
        ly_2 = gpuArray(l_2(2,:));
        lz_2 = gpuArray(l_2(3,:));
        dirvx_2 = gpuArray(dirv_2(1,:));
        dirvy_2 = gpuArray(dirv_2(2,:));
        dirvz_2 = gpuArray(dirv_2(3,:));
    end
end

if(dim == 2)
    if(sct_type == 3)
        ff_ms2D = parallel.gpu.CUDAKernel('ff_ms2D.ptx','ff_ms2D.cu');
        ff_ms2D.GridSize = [ceil((Nv * Nw)/ff_ms2D.MaxThreadsPerBlock) 1 1];
        ff_ms2D.ThreadBlockSize = [ff_ms2D.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_ms2D,'vSize',int32(Nv));
        setConstantMemory(ff_ms2D,'lSize',int32(Nw));
        setConstantMemory(ff_ms2D,'box_min',box_min);
        setConstantMemory(ff_ms2D,'box_max',box_max);
        setConstantMemory(ff_ms2D,'fw',ampfunc.forwardWeight);
        setConstantMemory(ff_ms2D,'sigt',sigt/2);
        setConstantMemory(ff_ms2D,'g_up',(1 - ampfunc.g * ampfunc.g)/(2*pi));
        setConstantMemory(ff_ms2D,'g_down1',1 + ampfunc.g * ampfunc.g);
        setConstantMemory(ff_ms2D,'g_down2',-2 * ampfunc.g);
    end
    
    if(sct_type == 2)
        ff_ms_tabulated2D = parallel.gpu.CUDAKernel('ff_ms_tabulated2D.ptx','ff_ms_tabulated2D.cu');
        ff_ms_tabulated2D.GridSize = [ceil((Nv * Nw)/ff_ms_tabulated2D.MaxThreadsPerBlock) 1 1];
        ff_ms_tabulated2D.ThreadBlockSize = [ff_ms_tabulated2D.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_ms_tabulated2D,'vSize',int32(Nv));
        setConstantMemory(ff_ms_tabulated2D,'lSize',int32(Nw));
        setConstantMemory(ff_ms_tabulated2D,'nAmpfunc',int32(numel(ampfunc.evalAmp)));
        setConstantMemory(ff_ms_tabulated2D,'box_min',box_min);
        setConstantMemory(ff_ms_tabulated2D,'box_max',box_max);
        setConstantMemory(ff_ms_tabulated2D,'sigt',sigt/2);
    end

    if(~noWl)
        ff_single2D = parallel.gpu.CUDAKernel('ff_single2D.ptx','ff_single2D.cu');
        ff_single2D.GridSize = [ceil((Nv * Nw)/ff_single2D.MaxThreadsPerBlock) 1 1];
        ff_single2D.ThreadBlockSize = [ff_single2D.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_single2D,'vSize',int32(Nv));
        setConstantMemory(ff_single2D,'lSize',int32(Nl));
        setConstantMemory(ff_single2D,'wSize',int32(Nw));
        setConstantMemory(ff_single2D,'box_min',box_min);
        setConstantMemory(ff_single2D,'box_max',box_max);
        setConstantMemory(ff_single2D,'sigt',sigt/2);
    else
        ff_single2D_noW = parallel.gpu.CUDAKernel('ff_single2D_noW.ptx','ff_single2D_noW.cu');
        ff_single2D_noW.GridSize = [ceil((Nv * Nl)/ff_single2D_noW.MaxThreadsPerBlock) 1 1];
        ff_single2D_noW.ThreadBlockSize = [ff_single2D_noW.MaxThreadsPerBlock 1 1];

        setConstantMemory(ff_single2D_noW,'vSize',int32(Nv));
        setConstantMemory(ff_single2D_noW,'lSize',int32(Nl));
        setConstantMemory(ff_single2D_noW,'box_min',box_min);
        setConstantMemory(ff_single2D_noW,'box_max',box_max);
        setConstantMemory(ff_single2D_noW,'sigt',sigt/2);
    end

    if(sct_type == 3)
        af_ang_vl_gpu2D = parallel.gpu.CUDAKernel('af_ang_vl_gpu2D.ptx','af_ang_vl_gpu2D.cu');
        af_ang_vl_gpu2D.GridSize = [ceil((Nv * Nl)/af_ang_vl_gpu2D.MaxThreadsPerBlock) 1 1];
        af_ang_vl_gpu2D.ThreadBlockSize = [af_ang_vl_gpu2D.MaxThreadsPerBlock 1 1];

        setConstantMemory(af_ang_vl_gpu2D,'vSize',int32(Nv));
        setConstantMemory(af_ang_vl_gpu2D,'lSize',int32(Nl));
        setConstantMemory(af_ang_vl_gpu2D,'fw',ampfunc.forwardWeight);
        setConstantMemory(af_ang_vl_gpu2D,'g_up',(1 - ampfunc.g * ampfunc.g)/(2*pi));
        setConstantMemory(af_ang_vl_gpu2D,'g_down1',1 + ampfunc.g * ampfunc.g);
        setConstantMemory(af_ang_vl_gpu2D,'g_down2',-2 * ampfunc.g);
    end
    
    if(sct_type == 2)
        af_ang_vl_gpu_tabulated2D = parallel.gpu.CUDAKernel('af_ang_vl_gpu_tabulated2D.ptx','af_ang_vl_gpu_tabulated2D.cu');
        af_ang_vl_gpu_tabulated2D.GridSize = [ceil((Nv * Nl)/af_ang_vl_gpu_tabulated2D.MaxThreadsPerBlock) 1 1];
        af_ang_vl_gpu_tabulated2D.ThreadBlockSize = [af_ang_vl_gpu_tabulated2D.MaxThreadsPerBlock 1 1];

        setConstantMemory(af_ang_vl_gpu_tabulated2D,'vSize',int32(Nv));
        setConstantMemory(af_ang_vl_gpu_tabulated2D,'lSize',int32(Nl));
        setConstantMemory(af_ang_vl_gpu_tabulated2D,'nAmpfunc',int32(numel(ampfunc.evalAmp)));
    end

    wl_1 = gpuArray(complex(wl_1));
    wv_1 = gpuArray(wv_1);
    vx_1 = gpuArray(v_1(1,:));
    vy_1 = gpuArray(v_1(2,:));
    lx_1 = gpuArray(l_1(1,:));
    ly_1 = gpuArray(l_1(2,:));
    dirvx_1 = gpuArray(dirv_1(1,:));
    dirvy_1 = gpuArray(dirv_1(2,:));

    if(correlationFlag)
        wl_2 = gpuArray(complex(wl_2));
        wv_2 = gpuArray(wv_2);
        vx_2 = gpuArray(v_2(1,:));
        vy_2 = gpuArray(v_2(2,:));
        lx_2 = gpuArray(l_2(1,:));
        ly_2 = gpuArray(l_2(2,:));
        dirvx_2 = gpuArray(dirv_2(1,:));
        dirvy_2 = gpuArray(dirv_2(2,:));
    end
    
end

if(sct_type == 2)
    if(dim == 2)
        evalAmp = gpuArray(complex(ampfunc.evalAmp));
    end
    
    if(dim == 3)
        evalAmp = gpuArray(ampfunc.evalAmp);
    end
end

%% Prepare for algorithm

u_0 = complex(zeros(Nv,Nw,'gpuArray'));

% Initiate output parameters
if(noWv)
    u = complex(zeros(Nv,Nw,'gpuArray'));
    us = complex(zeros(Nv,Nw,'gpuArray'));
else
    u = complex(zeros(size(wv_1,1),Nw,'gpuArray'));
    us = complex(zeros(size(wv_1,1),Nw,'gpuArray'));
end

% Box size
box_w = box_max-box_min;

% Pre-calculate single scattering rotation amplitude, only possible when
% both light and view are far field (otherwise it also dependent on the
% first scatter position)

if(dim == 3)
    af_ang_vl_1 = zeros(size(v_1,2), size(l_1,2), 'gpuArray');
    if(sct_type == 3)
        af_ang_vl_1 = feval(af_ang_vl_gpu, af_ang_vl_1, vx_1, vy_1, vz_1, lx_1, ly_1, lz_1);
    end
    
    if(sct_type == 2)
        af_ang_vl_1 = feval(af_ang_vl_gpu_tabulated, af_ang_vl_1, evalAmp, vx_1, vy_1, vz_1, lx_1, ly_1, lz_1);
    end
end

if(dim == 2)
    if(sct_type == 3)
        af_ang_vl_1 = zeros(size(v_1,2), size(l_1,2), 'gpuArray');
        af_ang_vl_1 = feval(af_ang_vl_gpu2D, af_ang_vl_1, vx_1, vy_1, lx_1, ly_1);
    end
    
    if(sct_type == 2)
        af_ang_vl_1 = gpuArray(complex(zeros(size(v_1,2), size(l_1,2))));
        af_ang_vl_1 = feval(af_ang_vl_gpu_tabulated2D, af_ang_vl_1, evalAmp, vx_1, vy_1, lx_1, ly_1);
    end
end

if(correlationFlag)
    if((mean(abs(v_1(:) - v_2(:))) + mean(abs(l_1(:) - l_2(:)))) < 1e-16 )
        af_ang_vl_2 = af_ang_vl_1;
    else
        if(dim == 3)
            af_ang_vl_2 = zeros(size(v_1,2), size(l_1,2), 'gpuArray');
            if(sct_type == 3)
                af_ang_vl_2 = feval(af_ang_vl_gpu, af_ang_vl_2, vx_2, vy_2, vz_2, lx_2, ly_2, lz_2);
            end

            if(sct_type == 2)
                af_ang_vl_2 = feval(af_ang_vl_gpu_tabulated, af_ang_vl_2, evalAmp, vx_2, vy_2, vz_2, lx_2, ly_2, lz_2);
            end
        end

        if(dim == 2)
            if(sct_type == 3)
                af_ang_vl_2 = zeros(size(v_1,2), size(l_1,2), 'gpuArray');
                af_ang_vl_2 = feval(af_ang_vl_gpu2D, af_ang_vl_2, vx_2, vy_2, lx_2, ly_2);
            end

            if(sct_type == 2)
                af_ang_vl_1 = gpuArray(complex(zeros(size(v_1,2), size(l_1,2))));
                af_ang_vl_2 = feval(af_ang_vl_gpu_tabulated2D, af_ang_vl_2, evalAmp, vx_2, vy_2, lx_2, ly_2);
            end
        end
    end

end

% threshold to begin kill particles with low weight
killThr=0.2;
V=prod(box_w);

% direction 0
d0 = zeros(dim,1); d0(end) = 1;

%% Begin the main loop
for itr=1:maxItr
%     if(mod(itr,1e3) == 0)
%         disp([num2str(round((itr/maxItr)*100)),'%'])
%     end
    
    if(zeroPercent ~= 0)
        curr_wl_1 = (rand(size(wl_1),'gpuArray') > zeroPercent) .* wl_1;
        curr_wv_1 = (rand(size(wv_1),'gpuArray') > zeroPercent) .* wv_1;
        
        if(correlationFlag)
            curr_wl_2 = (rand(size(wl_2),'gpuArray') > zeroPercent) .* wl_2;
            curr_wv_2 = (rand(size(wv_2),'gpuArray') > zeroPercent) .* wv_2;
        end
    else
        curr_wl_1 = wl_1;
        curr_wv_1 = wv_1;
        
        if(correlationFlag)
            curr_wl_2 = wl_2;
            curr_wv_2 = wv_2;
        end
    end
    
   %itr
    % Sample the first scatter
    % x: first scattering point
    % px: probability by which first point was sampled. Needed so that
    % importance sampling integration is weighted properly
    switch mod(smpFlg,10)
        case 1
            % uniform distribution
            x=rand(dim,1).*(box_w)+box_min; px=1;
        case 2
            % exponential distribution
            [x,px] = expSmpX(box_min,box_max,d0,sigt);
%             [x,px] = expSmpX([-10;-35],[10;35],d0,sigt(2));
        case 3
            error('bad!!')
        case 6
            x = ampfunc0.inX(:,1,1,itr);
            px = ampfunc0.inPxpw(1,1,itr);
    end
    
    % First scattering direction
    % w - sampled direction.
    % w0p - probability of the sampled direction, needed to compute inportance sampling integral correctly. 
    switch floor(smpFlg/10)
        case 1
            % uniform distribution
            w=randn(dim,1); w=w/norm(w);
            w0p=1/sqrt(2^(dim-1)*pi);
        case 2
            % g0
            w=smpampfunc_general(meanl, sct_type,ampfunc0);     
            w0p=(evalampfunc_general(meanl'*w,sct_type,ampfunc0,dim));
        case 3
            % multimode
            [w,w0p] = smpampfunc_multimode(l_1,sct_type,ampfunc0);
        case 6
            w = ampfunc0.inw(:,1,1,itr);
            w0p = ampfunc0.inPxpw(2,1,itr);
    end

    % rotation due to first scattering in multiple scattering case (s
    % function in article).
    af_l_1 = evalampfunc_general((w'*l_1),sct_type,ampfunc,dim)./w0p;
    
    % complex transmission (xi function in article) and attenuation term 
    % in multiple scattering (tau function in article) between first scattering
    % particle and the light source 
    e_l0_1 = evalphaseatt(x,l_1,sigt,lambda,box_min,box_max,1,dirl_1);
    e_l0_1 = gpuArray(e_l0_1);
    
    % complex volumetric throughput (ni function in article) of first
    % scattering event to be used with paths of length >1 (the multiple scattering
    % case)
    if(~noWl)
        e_l0_ms_1 = (sum(e_l0_1.*af_l_1.*curr_wl_1,2)).';
    else
        e_l0_ms_1 = (e_l0_1.*af_l_1).';
    end
    
    if(correlationFlag)
        % rotation due to first scattering in multiple scattering case (s
        % function in article).
        af_l_2 = evalampfunc_general((w'*l_2),sct_type,ampfunc,dim)./w0p;

        % complex transmission (xi function in article) and attenuation term 
        % in multiple scattering (tau function in article) between first scattering
        % particle and the light source 
        e_l0_2 = evalphaseatt(x,l_2,sigt,lambda,box_min,box_max,1,dirl_2);
        e_l0_2 = gpuArray(e_l0_2);

        % complex volumetric throughput (ni function in article) of first
        % scattering event to be used with paths of length >1 (the multiple scattering
        % case)
        if(~noWl)
            e_l0_ms_2 = (sum(e_l0_2.*af_l_2.*curr_wl_2,2)).';
        else
            e_l0_ms_2 = (e_l0_2.*af_l_2).';
        end
    end
      
    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
    weight=albedo;
    
    % begin paths sampling loop
    while 1

        pL=pL+1;

        % calculate the complex volumetric throughput for the last
        % scattering event in case of multiple scattering
        if (pL>1)
            constCont = sqrt(weight./px)*exp(2*pi*1i*rand);
    
            if(dim == 3)
                if(sct_type == 3)
                    setConstantMemory(ff_ms,'fastConstCopy',[x;ow;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_ms, u_0, e_l0_ms_1, vx_1, vy_1, vz_1, dirvx_1, dirvy_1, dirvz_1);
                end
                
                if(sct_type == 2)
                    setConstantMemory(ff_ms_tabulated,'fastConstCopy',[x;ow;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_ms_tabulated, u_0, e_l0_ms_1, evalAmp, vx_1, vy_1, vz_1, dirvx_1, dirvy_1, dirvz_1);
                end
            end
            
            if(dim == 2)
                if(sct_type == 3)
                    setConstantMemory(ff_ms2D,'fastConstCopy',[x;ow;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_ms2D, u_0, e_l0_ms_1, vx_1, vy_1, dirvx_1, dirvy_1);
                end
                
                if(sct_type == 2)
                    setConstantMemory(ff_ms_tabulated2D,'fastConstCopy',[x;ow;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_ms_tabulated2D, u_0, e_l0_ms_1, evalAmp, vx_1, vy_1, dirvx_1, dirvy_1);
                end
            end
            
            if(correlationFlag)
                if(dim == 3)
                    if(sct_type == 3)
                        u_2 = curr_wv_2 * feval(ff_ms, u_0, e_l0_ms_2, vx_2, vy_2, vz_2, dirvx_2, dirvy_2, dirvz_2);
                    end
                    
                    if(sct_type == 2)
                        u_2 = curr_wv_2 * feval(ff_ms_tabulated, u_0, e_l0_ms_2, evalAmp, vx_2, vy_2, vz_2, dirvx_2, dirvy_2, dirvz_2);
                    end
                end
                
                if(dim == 2)
                    if(sct_type == 3)
                        u_2 = curr_wv_2 * feval(ff_ms2D, u_0, e_l0_ms_2, vx_2, vy_2, dirvx_2, dirvy_2);
                    end
                    
                    if(sct_type == 2)
                        u_2 = curr_wv_2 * feval(ff_ms_tabulated2D, u_0, e_l0_ms_2, evalAmp, vx_2, vy_2, dirvx_2, dirvy_2);
                    end
                end
                
                u = u + u_1 .* conj(u_2);
            else
                u = u + u_1;
            end
        end

        % Update field with next-event estimation
        if (pL==1)
            constCont = sqrt(weight./px)*exp(2*pi*1i*rand);
            
            if(dim == 3)
                if(~noWl)
                    setConstantMemory(ff_single,'fastConstCopy',[x;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_single, u_0, curr_wl_1, e_l0_1, af_ang_vl_1, ...
                        vx_1, vy_1, vz_1, dirvx_1, dirvy_1, dirvz_1);
                else
                    setConstantMemory(ff_single_noW,'fastConstCopy',[x;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_single_noW, u_0, e_l0_1, af_ang_vl_1, ...
                        vx_1, vy_1, vz_1, dirvx_1, dirvy_1, dirvz_1);
                end
            end
            
            if(dim == 2)
                if(~noWl)
                    setConstantMemory(ff_single2D,'fastConstCopy',[x;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_single2D, u_0, curr_wl_1, e_l0_1, real(af_ang_vl_1), ...
                        vx_1, vy_1, dirvx_1, dirvy_1);
                else
                    setConstantMemory(ff_single2D_noW,'fastConstCopy',[x;real(constCont);imag(constCont)]);
                    u_1 = curr_wv_1 * feval(ff_single2D_noW, u_0, e_l0_1, af_ang_vl_1, ...
                        vx_1, vy_1, dirvx_1, dirvy_1);
                end
            end
            
            if(correlationFlag)
                if(dim == 3)
                    u_2 = curr_wv_2 * feval(ff_single, u_0, curr_wl_2, e_l0_2, af_ang_vl_2, ...
                        vx_2, vy_2, vz_2, dirvx_2, dirvy_2, dirvz_2);
                end
                
                if(dim == 2)
                    if(~noWl)
                        u_2 = curr_wv_2 * feval(ff_single2D, u_0, curr_wl_2, e_l0_2, real(af_ang_vl_2), ...
                            vx_2, vy_2, dirvx_2, dirvy_2);
                    else
                        u_2 = curr_wv_2 * feval(ff_single2D_noW, u_0, e_l0_2, af_ang_vl_2, ...
                            vx_2, vy_2, dirvx_2, dirvy_2);
                    end
                end
                
                us = us + u_1 .* conj(u_2);
                u = u + u_1 .* conj(u_2);
            else
                us = us + u_1;
                u = u + u_1;
            end
        end

        % advance to the next scattering event
        if(smpFlg == 66)
            x = ampfunc0.inX(:,pL+1,1,itr);
        else
            d=-log(-rand+1)/(sigt);
            x=x+d*w;
        end

        % move to the next particle if the next scattering event is outside
        % the box
        if(max(x>box_max) || max(x<box_min))
            break
        end

        % albedo intensity reduction. If the intensity is too low, it is
        % possible that the particle will be killed
        if weight<killThr
            if rand>albedo
                break
            end
        else
            weight=weight*albedo;
        end

        % Sample new scatteing direction
        ow=w;
        
        if(smpFlg == 66)
            w = ampfunc0.inw(:,pL+1,1,itr);
        else
            w=smpampfunc_general(ow, sct_type,ampfunc);
        end

    end
  
end

%% Normalization

if(~correlationFlag)
    u=u*sqrt(1/maxItr*V*sigt);
    us=us*sqrt(1/maxItr*V*sigt);
else
    u=u*(1/maxItr*V*sigt);
    us=us*(1/maxItr*V*sigt);
    
%     u=u*(1/maxItr);
%     us=us*(1/maxItr);
end

u = gather(u);
us = gather(us);

end
