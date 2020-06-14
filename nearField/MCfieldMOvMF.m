function [u,us,u_nomean,xRep,wRep,xVec,wVec,pxpwVec] = MCfieldMOvMF(sigt,albedo,box_min,box_max,maxItr,pathsNum,lambda,smpFlg,smpFunc,movmf,sct_type,ampfunc,gpuEnable, ...
    apertureVmf_l,apertureVmf_v)
u_nomean = 0;
xRep = 0;
wRep = 0;

outxpos = (nargout > 5);
maxBounces = 100;

if(outxpos)
    xVec = inf(3,maxBounces,pathsNum,maxItr);
    wVec = inf(3,maxBounces,pathsNum,maxItr);
    pxpwVec = inf(2,pathsNum,maxItr);
end

dim = size(box_min,1);

u_size = max(apertureVmf_l.dim,apertureVmf_v.dim);
u_size = u_size(2:end);
u_size(end+1) = pathsNum;

u_output_size = max(apertureVmf_l.dim,apertureVmf_v.dim);
u_output_size = u_output_size(2:4);

dirl = apertureVmf_l.dir;
dirv = apertureVmf_v.dir;

correlationFlag = (u_size(4) == 2);
if(correlationFlag) 
    u_size(4) = 1;
end


%% Gpu function
if(gpuEnable)
    lDim = apertureVmf_l.dim;
    lDim(1) = movmf.dim(1);
    lDim(apertureVmf_v.dirDim) = apertureVmf_v.dim(apertureVmf_v.dirDim);

    if(correlationFlag)
        gpuFunc.integrateMult = parallel.gpu.CUDAKernel('integrateMult_correlation.ptx','integrateMult_correlation.cu');
    else
        gpuFunc.integrateMult = parallel.gpu.CUDAKernel('integrateMult_render.ptx','integrateMult_render.cu');
    end
    
    gpuFunc.integrateMult.GridSize = [ceil(prod(u_size)/gpuFunc.integrateMult.MaxThreadsPerBlock) 1 1];
    gpuFunc.integrateMult.ThreadBlockSize = [gpuFunc.integrateMult.MaxThreadsPerBlock 1 1];

    setConstantMemory(gpuFunc.integrateMult,'box_min',box_min);
    setConstantMemory(gpuFunc.integrateMult,'box_max',box_max);
    setConstantMemory(gpuFunc.integrateMult,'sigt',sigt/2);

    setConstantMemory(gpuFunc.integrateMult,'dirDimNum',int32(apertureVmf_v.dirDim));
    setConstantMemory(gpuFunc.integrateMult,'uDimProd',int32(cumprod(u_size)));
    setConstantMemory(gpuFunc.integrateMult,'lDim',int32(lDim));
    setConstantMemory(gpuFunc.integrateMult,'vDim',int32(apertureVmf_v.dim));
    setConstantMemory(gpuFunc.integrateMult,'lMixtureAlpha',movmf.alpha);

    if(correlationFlag)
        gpuFunc.MSintegrateMult = parallel.gpu.CUDAKernel('MSintegrateMult_correlation.ptx','MSintegrateMult_correlation.cu');
    else
        gpuFunc.MSintegrateMult = parallel.gpu.CUDAKernel('MSintegrateMult_render.ptx','MSintegrateMult_render.cu');
    end
    
    gpuFunc.MSintegrateMult.GridSize = [ceil(prod(u_size)/gpuFunc.MSintegrateMult.MaxThreadsPerBlock) 1 1];
    gpuFunc.MSintegrateMult.ThreadBlockSize = [gpuFunc.MSintegrateMult.MaxThreadsPerBlock 1 1];

    setConstantMemory(gpuFunc.MSintegrateMult,'box_min',box_min);
    setConstantMemory(gpuFunc.MSintegrateMult,'box_max',box_max);
    setConstantMemory(gpuFunc.MSintegrateMult,'sigt',sigt/2);

    setConstantMemory(gpuFunc.MSintegrateMult,'dirDimNum',int32(apertureVmf_v.dirDim));
    setConstantMemory(gpuFunc.MSintegrateMult,'uDimProd',int32(cumprod(u_size)));
    setConstantMemory(gpuFunc.MSintegrateMult,'lDim',int32(apertureVmf_l.dim));
    setConstantMemory(gpuFunc.MSintegrateMult,'vDim',int32(apertureVmf_v.dim));
    setConstantMemory(gpuFunc.MSintegrateMult,'lMixtureAlpha',movmf.alpha);
    setConstantMemory(gpuFunc.MSintegrateMult,'mixturesNum',int32(movmf.dim(1)));
    
    movmf_gpu = movmfToGpu(movmf);
    apertureVmf_v_gpu = movmfToGpu(apertureVmf_v);
    dirv_gpu = gpuArray(dirv);

end
%% Prepare for algorithm
if(gpuEnable)
    u = complex(zeros(u_output_size,'gpuArray'));
    us = complex(zeros(u_output_size,'gpuArray'));
else
    u = complex(zeros(u_output_size));
    us = complex(zeros(u_output_size));
end

% Box size
box_w = box_max-box_min;
V=prod(box_w);

% threshold to begin kill particles with low weight
killThr=0.2;

%% Begin the main loop
for itr=1:maxItr
%     if(mod(itr,1e1) == 0) 
%         disp([num2str(round(1e2 * itr/maxItr)), '%'])
%     end
    % Sample the first scatter
    % x: first scattering point
    % px: probability by which first point was sampled. Needed so that
    % importance sampling integration is weighted properly
    switch mod(smpFlg,10)
        case 1
            % uniform distribution
            x=rand(dim,1,1,1,1,pathsNum).*(box_w)+box_min;
            px=ones(1,1,1,1,1,pathsNum);
        case 2
            % exponential distribution
            [x,px] = expSmpX(box_min,box_max,[0;0;1],sigt,pathsNum);
        case 3
%             [x,px,xMissIter] = smpVmfBeam(apertureVmf_l,smpFunc,box_min,box_max);
%             px = px * V; % The reson I multiple with V - same prob as uniform
%             xRep = xRep + xMissIter;
        case 4
            [x,px,xMissIter,n] = smpVmfBeamSum(apertureVmf_l,smpFunc,box_min,box_max,pathsNum);
            px = px * V; % The reson I multiple with V - same prob as uniform
            xRep = xRep + xMissIter;
            
        case 7
            x = smpFunc;
            px = 1;
    end
%     px = max(px,smpFunc(1).pxMin);
    
    dz = cubeDist(x,box_min,box_max,-dirl);
    throughputVmf_l = movmfThroughput(apertureVmf_l,x,-1,sigt,dz);
    conv_vmf_l0 = movmfConv3(throughputVmf_l,movmf,dirv,apertureVmf_v.dirDim);
    
    % First scattering direction
    %implement an option to sample w isotropically (add apropriate input flag).
    %otherwise implement an option that samples one of the incoming
    %illumination of interest and sample a direction from that Gaussian. 
    %sampling probability w0p should be
    switch floor(smpFlg/10)
        case 1
            w=randn(dim,1,1,1,1,pathsNum); w=w./sqrt(sum(w.^2,1));
            w0p=1./sqrt(2^(dim-1)*pi);
        case 3
%             [w,w0p] = vmfFirstDirection(conv_vmf_l0);
        case 4
%             n = randi(prod(conv_vmf_l0.dim(2:end)),[1,1,1,1,pathsNum]);
%             [w,w0p] = vmfFirstDirectionSum_sameBeam(conv_vmf_l0,n);
        case 5
            [w,w0p] = vmfFirstDirectionSum_sameBeam(conv_vmf_l0,n);
    end
    
    if(outxpos)
        pxpwVec(1,:,itr) = px;
        pxpwVec(2,:,itr) = w0p;
    end
    
    % for each lighting direction, the viewing direction gets different
    % direction
    
    if(gpuEnable)
%         randPhase = exp(2*pi*1i*rand(1,1,1,1,1,pathsNum,'gpuArray'));
        randPhase = gpuArray(exp(2*pi*1i*rand(1,1,1,1,1,pathsNum)));
        conv_vmf_l0_gpu = movmfToGpu(conv_vmf_l0);
        
        us = feval(gpuFunc.integrateMult,us,dirv_gpu,gpuArray(x),(2*pi*randPhase./sqrt(px)), ...
            conv_vmf_l0_gpu.mu1, conv_vmf_l0_gpu.mu2, conv_vmf_l0_gpu.mu3, conv_vmf_l0_gpu.c, ...
            apertureVmf_v_gpu.mu1, apertureVmf_v_gpu.mu2, apertureVmf_v_gpu.mu3, apertureVmf_v_gpu.c);
    else
        randPhase = exp(2*pi*1i*rand(1,1,1,1,1,pathsNum));
        dz = cubeDist(x,box_min,box_max,dirv);
        throughputVmf_v = movmfThroughput(apertureVmf_v,x,1,sigt,dz);
        movmf_mult = movmfMultiple(conv_vmf_l0,throughputVmf_v,false);
        [integrateMult] = movmfIntegrate(movmf_mult);
        us = us + squeeze(integrateMult .* randPhase / sqrt(px));
    end
    
    rotatedMovmf_l = movmfRotate(movmf,w);
    
    throughputVmf_l_times_movmf = movmfMultiple(throughputVmf_l,rotatedMovmf_l,true);
    throughputVmf_l_times_movmf.c = throughputVmf_l_times_movmf.c - log(w0p);
    el = movmfIntegrate(throughputVmf_l_times_movmf);

    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
%     weight=albedo;
    
    if(gpuEnable)
        activatedPaths = true(pathsNum,1,'gpuArray');
        el_gpu = gpuArray(el);
    end
    
    % begin paths sampling loop
    while pL < maxBounces

        pL=pL+1;
        
        if(outxpos)
            xVec(:,pL,:,itr) = x;
            wVec(:,pL,:,itr) = w;
        end


        % calculate the complex volumetric throughput for the last
        % scattering event in case of multiple scattering
        if (pL>1)
                
            if(gpuEnable)
%                 randPhase = exp(2*pi*1i*rand(1,1,1,1,1,pathsNum,'gpuArray'));
                randPhase = gpuArray(exp(2*pi*1i*rand(1,1,1,1,1,pathsNum)));

                % rotate the scattering function on direction of ow
                rotatedMovmf_v = movmfRotate(movmf_gpu,ow);
                
                u = feval(gpuFunc.MSintegrateMult,u, el_gpu, dirv, ...
                    gpuArray(x), (2 * pi * randPhase ./ sqrt(px)), activatedPaths, ...
                    apertureVmf_v.mu1, apertureVmf_v.mu2, apertureVmf_v.mu3, apertureVmf_v.c, ...
                    rotatedMovmf_v.mu1, rotatedMovmf_v.mu2, rotatedMovmf_v.mu3, rotatedMovmf_v.c);
            else
                randPhase = exp(2*pi*1i*rand(1,1,1,1,1,pathsNum));

                % rotate the scattering function on direction of ow
                rotatedMovmf_v = movmfRotate(movmf,ow);
                
                dz = cubeDist(x,box_min,box_max,dirv);
                throughputVmf_v = movmfThroughput(apertureVmf_v,x,1,sigt,dz);
                throughputVmf_v_times_movmf = movmfMultiple(throughputVmf_v,rotatedMovmf_v,true);
                ev = movmfIntegrate(throughputVmf_v_times_movmf);
                u_res = squeeze(ev .* el ./ sqrt(px) .* randPhase);
                u(:,:,:,activatedPaths) = u(:,:,:,activatedPaths) + u_res(:,:,:,activatedPaths);
            end
        end

       
        % advance to the next scattering event
        d=-log(-rand(1,1,1,1,1,pathsNum)+1)/(sigt);
        x = x+d.*w;

        % move to the next particle if the next scattering event is outside
        % the box
        outsideIdx = ...
            x(1,1,1,1,1,:) > box_max(1) | x(1,1,1,1,1,:) < box_min(1) | ...
            x(2,1,1,1,1,:) > box_max(2) | x(2,1,1,1,1,:) < box_min(2) | ...
            x(3,1,1,1,1,:) > box_max(3) | x(3,1,1,1,1,:) < box_min(3) ;
        
        activatedPaths(outsideIdx) = false;
        
        if(~any(activatedPaths))
            break
        end

%         % albedo intensity reduction. If the intensity is too low, it is
%         % possible that the particle will be killed
%         if weight<killThr
%             if rand>albedo
%                 break
%             end
%         else
%             weight=weight*albedo;
%         end

        % Sample new scatteing direction
        ow=w;
        
        w = smpampfunc_general(ow, sct_type,ampfunc,pathsNum);
    end
  
end

%% Normalization
% u=u*sqrt(1/maxItr*V*sigt(2));
% us=us*sqrt(1/maxItr*V*sigt(2));

if(~correlationFlag)
    u=u*sqrt(V*sigt / (maxItr*pathsNum + xRep));
    us=us*sqrt(V*sigt / (maxItr*pathsNum + xRep));
else
    u=u*(V*sigt / (maxItr*pathsNum + xRep));
    us=us*(V*sigt / (maxItr*pathsNum + xRep));
end

u = u+us;

u = gather(u);
us = gather(us);

end
