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
correlationFlag = (size(apertureVmf_l.mu1,2) == 2);

%% Calculate dims
MAX_DIM = 32;

conv_maxDim = max( ...
    [ndims(apertureVmf_l.mu1) , ndims(apertureVmf_l.mu2) , ndims(apertureVmf_l.mu3) , ...
     ndims(apertureVmf_v.dir1), ndims(apertureVmf_v.dir2), ndims(apertureVmf_v.dir3), ...
     ndims(lambda)]);
 
% convolution dims
conv_mu1_size = size(apertureVmf_l.mu1); conv_mu1_size(end+1:conv_maxDim) = 1;
conv_mu2_size = size(apertureVmf_l.mu2); conv_mu2_size(end+1:conv_maxDim) = 1;
conv_mu3_size = size(apertureVmf_l.mu3); conv_mu3_size(end+1:conv_maxDim) = 1;
conv_c_size   = size(apertureVmf_l.c)  ; conv_c_size(end+1:conv_maxDim)   = 1;

conv_th_mu1_size = conv_mu1_size; conv_th_mu1_size(3) = pathsNum;
conv_th_mu2_size = conv_mu2_size; conv_th_mu2_size(3) = pathsNum;
conv_th_mu3_size = conv_mu3_size; conv_th_mu3_size(3) = pathsNum;
conv_th_c_size   = conv_c_size  ; conv_th_c_size(3)   = pathsNum;

conv_vdir1_size = size(apertureVmf_v.dir1); conv_vdir1_size(end+1:conv_maxDim) = 1;
conv_vdir2_size = size(apertureVmf_v.dir2); conv_vdir2_size(end+1:conv_maxDim) = 1;
conv_vdir3_size = size(apertureVmf_v.dir3); conv_vdir3_size(end+1:conv_maxDim) = 1;

conv_ldir1_size = size(apertureVmf_l.dir1); conv_ldir1_size(end+1:conv_maxDim) = 1;
conv_ldir2_size = size(apertureVmf_l.dir2); conv_ldir2_size(end+1:conv_maxDim) = 1;
conv_ldir3_size = size(apertureVmf_l.dir3); conv_ldir3_size(end+1:conv_maxDim) = 1;

conv_lambda_size = size(lambda); conv_lambda_size(end+1:conv_maxDim) = 1;

convDim = max([ ...
    conv_mu1_size   ; conv_mu2_size   ; conv_mu3_size    ; ...
    conv_vdir1_size ; conv_vdir2_size ; conv_vdir3_size  ; ...
    conv_lambda_size]);

convDim(3) = pathsNum;
convDim(1) = numel(movmf.alpha);

convDimProd = zeros(15,MAX_DIM);
convDimProd(1 ,1:numel(conv_mu1_size))    = cumprod(conv_mu1_size)   ;
convDimProd(2 ,1:numel(conv_mu2_size))    = cumprod(conv_mu2_size)   ;
convDimProd(3 ,1:numel(conv_mu3_size))    = cumprod(conv_mu3_size)   ;
convDimProd(4 ,1:numel(conv_c_size))      = cumprod(conv_c_size)     ;
convDimProd(5 ,1:numel(conv_vdir1_size))  = cumprod(conv_vdir1_size) ;
convDimProd(6 ,1:numel(conv_vdir2_size))  = cumprod(conv_vdir2_size) ;
convDimProd(7 ,1:numel(conv_vdir3_size))  = cumprod(conv_vdir3_size) ;
convDimProd(8 ,1:numel(conv_th_mu1_size)) = cumprod(conv_th_mu1_size);
convDimProd(9 ,1:numel(conv_th_mu2_size)) = cumprod(conv_th_mu2_size);
convDimProd(10,1:numel(conv_th_mu3_size)) = cumprod(conv_th_mu3_size);
convDimProd(11,1:numel(conv_th_c_size))   = cumprod(conv_th_c_size)  ;
convDimProd(12,1:numel(conv_ldir1_size))  = cumprod(conv_ldir1_size) ;
convDimProd(13,1:numel(conv_ldir2_size))  = cumprod(conv_ldir2_size) ;
convDimProd(14,1:numel(conv_ldir3_size))  = cumprod(conv_ldir3_size) ;
convDimProd(15,1:numel(conv_lambda_size)) = cumprod(conv_lambda_size) ;

zerosIdx = ([ones(15,1),diff(convDimProd,[],2)] ~= 0);

convDimProd = circshift(convDimProd,1,2);
convDimProd = zerosIdx .* convDimProd;

% scattering dims
scatt_maxDim = max( ...
    [ndims(apertureVmf_v.mu1) , ndims(apertureVmf_v.mu2) , ndims(apertureVmf_v.mu3) , ...
     numel(convDim), ndims(lambda)]);
 
scatt_l_size = convDim; scatt_l_size(end+1:scatt_maxDim) = 1;

scatt_mu1_size = size(apertureVmf_v.mu1); scatt_mu1_size(end+1:scatt_maxDim) = 1;
scatt_mu2_size = size(apertureVmf_v.mu2); scatt_mu2_size(end+1:scatt_maxDim) = 1;
scatt_mu3_size = size(apertureVmf_v.mu3); scatt_mu3_size(end+1:scatt_maxDim) = 1;
scatt_c_size   = size(apertureVmf_v.c)  ; scatt_c_size(end+1:scatt_maxDim)   = 1;

scatt_vdir1_size = size(apertureVmf_v.dir1); scatt_vdir1_size(end+1:scatt_maxDim) = 1;
scatt_vdir2_size = size(apertureVmf_v.dir2); scatt_vdir2_size(end+1:scatt_maxDim) = 1;
scatt_vdir3_size = size(apertureVmf_v.dir3); scatt_vdir3_size(end+1:scatt_maxDim) = 1;

scatt_lambda_size = size(lambda); scatt_lambda_size(end+1:scatt_maxDim) = 1;

u_size = max([ ...
    scatt_l_size
    scatt_mu1_size   ; scatt_mu2_size   ; scatt_mu3_size    ; ...
    scatt_vdir1_size ; scatt_vdir2_size ; scatt_vdir3_size  ; ...
    scatt_lambda_size]);

scattThreadsDims = u_size(3:end);
u_size = u_size(4:end);

scattDimProd = zeros(9,MAX_DIM);
scattDimProd(1,1:numel(scatt_l_size))      = cumprod(scatt_l_size)     ;
scattDimProd(2,1:numel(scatt_mu1_size))    = cumprod(scatt_mu1_size)   ;
scattDimProd(3,1:numel(scatt_mu2_size))    = cumprod(scatt_mu2_size)   ;
scattDimProd(4,1:numel(scatt_mu3_size))    = cumprod(scatt_mu3_size)   ;
scattDimProd(5,1:numel(scatt_c_size))      = cumprod(scatt_c_size)     ;
scattDimProd(6,1:numel(scatt_vdir1_size))  = cumprod(scatt_vdir1_size) ;
scattDimProd(7,1:numel(scatt_vdir2_size))  = cumprod(scatt_vdir2_size) ;
scattDimProd(8,1:numel(scatt_vdir3_size))  = cumprod(scatt_vdir3_size) ;
scattDimProd(9,1:numel(scatt_lambda_size)) = cumprod(scatt_lambda_size);

zerosIdx = ([ones(9,1),diff(scattDimProd,[],2)] ~= 0);
scattDimProd = circshift(scattDimProd,1,2);
scattDimProd = zerosIdx .* scattDimProd;
scattDimProd = circshift(scattDimProd,-2,2);

% multiple scattering dims
mscatt_l_size = [1,scatt_l_size(2:end)];

mscattDimProd = zeros(9,MAX_DIM);
mscattDimProd(1,1:numel(mscatt_l_size))     = cumprod(mscatt_l_size)    ;
mscattDimProd(2,1:numel(scatt_mu1_size))    = cumprod(scatt_mu1_size)   ;
mscattDimProd(3,1:numel(scatt_mu2_size))    = cumprod(scatt_mu2_size)   ;
mscattDimProd(4,1:numel(scatt_mu3_size))    = cumprod(scatt_mu3_size)   ;
mscattDimProd(5,1:numel(scatt_c_size))      = cumprod(scatt_c_size)     ;
mscattDimProd(6,1:numel(scatt_vdir1_size))  = cumprod(scatt_vdir1_size) ;
mscattDimProd(7,1:numel(scatt_vdir2_size))  = cumprod(scatt_vdir2_size) ;
mscattDimProd(8,1:numel(scatt_vdir3_size))  = cumprod(scatt_vdir3_size) ;
mscattDimProd(9,1:numel(scatt_lambda_size)) = cumprod(scatt_lambda_size);

zerosIdx = ([ones(9,1),diff(mscattDimProd,[],2)] ~= 0);
mscattDimProd = circshift(mscattDimProd,1,2);
mscattDimProd = zerosIdx .* mscattDimProd;
mscattDimProd = circshift(mscattDimProd,-2,2);

% el dims
elDim = convDim; elDim(1) = 1;

elDimProd = zeros(4,MAX_DIM);
elDimProd(1,1:numel(conv_th_mu1_size)) = cumprod(conv_th_mu1_size);
elDimProd(2,1:numel(conv_th_mu2_size)) = cumprod(conv_th_mu2_size);
elDimProd(3,1:numel(conv_th_mu3_size)) = cumprod(conv_th_mu3_size);
elDimProd(4,1:numel(conv_th_c_size))   = cumprod(conv_th_c_size)  ;

zerosIdx = ([ones(4,1),diff(elDimProd,[],2)] ~= 0);

elDimProd = circshift(elDimProd,1,2);
elDimProd = zerosIdx .* elDimProd;


%% Gpu function
if(gpuEnable)
    gpuFunc.movmfConv = parallel.gpu.CUDAKernel('nf_core.ptx','nf_core.cu','movmfConv');
    gpuFunc.movmfConv.GridSize = [ceil(prod(convDim)/gpuFunc.movmfConv.MaxThreadsPerBlock) 1 1];
    gpuFunc.movmfConv.ThreadBlockSize = [gpuFunc.movmfConv.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.movmfConv,'mixtureMu',movmf.mu3);
    setConstantMemory(gpuFunc.movmfConv,'mixtureC',movmf.c);
    setConstantMemory(gpuFunc.movmfConv,'conv_dimProd',uint32(cumprod(convDim)));
    setConstantMemory(gpuFunc.movmfConv,'conv_maxDim',uint32(conv_maxDim));
    setConstantMemory(gpuFunc.movmfConv,'conv_inDimProd',uint32(convDimProd.'));
    setConstantMemory(gpuFunc.movmfConv,'box_min',box_min);
    setConstantMemory(gpuFunc.movmfConv,'box_max',box_max);
    setConstantMemory(gpuFunc.movmfConv,'sigt',sigt/2    );
   
    gpuFunc.singleScattering = parallel.gpu.CUDAKernel('nf_core.ptx','nf_core.cu','singleScattering');
    gpuFunc.singleScattering.GridSize = [ceil(prod(scattThreadsDims)/gpuFunc.singleScattering.MaxThreadsPerBlock) 1 1];
    gpuFunc.singleScattering.ThreadBlockSize = [gpuFunc.singleScattering.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.singleScattering,'mixtureAlpha',movmf.alpha);
    setConstantMemory(gpuFunc.singleScattering,'mixturesNum',uint32(numel(movmf.alpha)));
    setConstantMemory(gpuFunc.singleScattering,'scatt_dimProd',uint32(cumprod(scattThreadsDims)));
    setConstantMemory(gpuFunc.singleScattering,'scatt_maxDim',uint32(scatt_maxDim - 2));
    setConstantMemory(gpuFunc.singleScattering,'scatt_inDimProd',uint32(scattDimProd.'));
    setConstantMemory(gpuFunc.singleScattering,'corrParamNum',uint32(scatt_c_size(2)));
    
    gpuFunc.multipleScattering = parallel.gpu.CUDAKernel('nf_core.ptx','nf_core.cu','multipleScattering');
    gpuFunc.multipleScattering.GridSize = [ceil(prod(scattThreadsDims)/gpuFunc.multipleScattering.MaxThreadsPerBlock) 1 1];
    gpuFunc.multipleScattering.ThreadBlockSize = [gpuFunc.multipleScattering.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.multipleScattering,'mscatt_inDimProd',uint32(mscattDimProd.'));
    
    gpuFunc.integrateEl = parallel.gpu.CUDAKernel('nf_core.ptx','nf_core.cu','integrateEl');
    gpuFunc.integrateEl.GridSize = [ceil(prod(elDim)/gpuFunc.integrateEl.MaxThreadsPerBlock) 1 1];
    gpuFunc.integrateEl.ThreadBlockSize = [gpuFunc.integrateEl.MaxThreadsPerBlock 1 1];
    
    setConstantMemory(gpuFunc.integrateEl,'el_dimProd',uint32(cumprod(elDim)));
    setConstantMemory(gpuFunc.integrateEl,'el_inDimProd',uint32(elDimProd.'));
    
    conv_vmf_l0.mu1 = gpuArray(complex(zeros(convDim)));
    conv_vmf_l0.mu2 = gpuArray(complex(zeros(convDim)));
    conv_vmf_l0.mu3 = gpuArray(complex(zeros(convDim)));
    conv_vmf_l0.c = gpuArray(complex(zeros(convDim)));
    
    apertureVmf_v.mu1 = gpuArray(complex(apertureVmf_v.mu1));
    apertureVmf_v.mu2 = gpuArray(complex(apertureVmf_v.mu2));
    apertureVmf_v.mu3 = gpuArray(complex(apertureVmf_v.mu3));
    apertureVmf_v.c = gpuArray(complex(apertureVmf_v.c));
    
    apertureVmf_v.dir1 = gpuArray(apertureVmf_v.dir1);
    apertureVmf_v.dir2 = gpuArray(apertureVmf_v.dir2);
    apertureVmf_v.dir3 = gpuArray(apertureVmf_v.dir3);
    
    apertureVmf_l.mu1 = gpuArray(complex(apertureVmf_l.mu1));
    apertureVmf_l.mu2 = gpuArray(complex(apertureVmf_l.mu2));
    apertureVmf_l.mu3 = gpuArray(complex(apertureVmf_l.mu3));
    apertureVmf_l.c = gpuArray(complex(apertureVmf_l.c));
    
    apertureVmf_l.dir1 = gpuArray(apertureVmf_l.dir1);
    apertureVmf_l.dir2 = gpuArray(apertureVmf_l.dir2);
    apertureVmf_l.dir3 = gpuArray(apertureVmf_l.dir3);
    
    throughputVmf_l.mu1 = gpuArray(complex(zeros(conv_th_mu1_size)));
    throughputVmf_l.mu2 = gpuArray(complex(zeros(conv_th_mu2_size)));
    throughputVmf_l.mu3 = gpuArray(complex(zeros(conv_th_mu3_size)));
    throughputVmf_l.c   = gpuArray(complex(zeros(conv_th_c_size)))  ;
    
    el_zeros = gpuArray(complex(zeros(elDim)));
    lambda = gpuArray(lambda);
end
%% Prepare for algorithm
if(gpuEnable)
    u = complex(zeros(u_size,'gpuArray'));    
    us = complex(zeros(u_size,'gpuArray'));    
else
    u = complex(zeros(u_size));
    us = complex(zeros(u_size));
end

% Box size
box_w = box_max-box_min;
V=prod(box_w);

% threshold to begin kill particles with low weight
killThr=0.2;

dirl_zeros = 0 * apertureVmf_l.dir1 + 0 * apertureVmf_l.dir2 + 0 * apertureVmf_l.dir3;
dirl = cat(1,apertureVmf_l.dir1 + dirl_zeros,apertureVmf_l.dir2 + dirl_zeros,apertureVmf_l.dir3 + dirl_zeros);

if(~gpuEnable)
    dirv_zeros = 0 * apertureVmf_v.dir1 + 0 * apertureVmf_v.dir2 + 0 * apertureVmf_v.dir3;
    dirv = cat(1,apertureVmf_v.dir1 + dirv_zeros,apertureVmf_v.dir2 + dirv_zeros,apertureVmf_v.dir3 + dirv_zeros);
end

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
            x=rand(dim,1,pathsNum).*(box_w)+box_min;
            px=ones(1,1,pathsNum);
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
    
    if(gpuEnable)           
        [conv_vmf_l0.mu1,conv_vmf_l0.mu2,conv_vmf_l0.mu3,conv_vmf_l0.c, ...
            throughputVmf_l.mu1,throughputVmf_l.mu2,throughputVmf_l.mu3,throughputVmf_l.c] = ...
            feval(gpuFunc.movmfConv,conv_vmf_l0.mu1,conv_vmf_l0.mu2,conv_vmf_l0.mu3,conv_vmf_l0.c,...
            throughputVmf_l.mu1,throughputVmf_l.mu2,throughputVmf_l.mu3,throughputVmf_l.c, lambda, x, ...
            apertureVmf_l.mu1, apertureVmf_l.mu2, apertureVmf_l.mu3, apertureVmf_l.c, ...
            apertureVmf_v.dir1,apertureVmf_v.dir2,apertureVmf_v.dir3, ...
            apertureVmf_l.dir1,apertureVmf_l.dir2,apertureVmf_l.dir3);
    else
        dz = cubeDist(x,box_min,box_max,-dirl);
        throughputVmf_l = movmfThroughput(apertureVmf_l,lambda,x,-1,sigt,dz);
        
        conv_vmf_l0 = movmfConv(throughputVmf_l,apertureVmf_v,movmf);
    end
    
    % First scattering direction
    %implement an option to sample w isotropically (add apropriate input flag).
    %otherwise implement an option that samples one of the incoming
    %illumination of interest and sample a direction from that Gaussian. 
    %sampling probability w0p should be
    switch floor(smpFlg/10)
        case 1
            w=randn(dim,1,pathsNum); w=w./sqrt(sum(w.^2,1));
            w0p=ones(1,1,pathsNum)./sqrt(2^(dim-1)*pi);
        case 3
%             [w,w0p] = vmfFirstDirection(conv_vmf_l0);
        case 4
            Naperture = numel(apertureVmf_l.c);
            n = randi(Naperture,[1,1,pathsNum]);
            [w,w0p] = vmfFirstDirectionSum_sameBeam(conv_vmf_l0,movmf.alpha,n);
        case 5
            [w,w0p] = vmfFirstDirectionSum_sameBeam(conv_vmf_l0,movmf.alpha,n);
    end
    
%     if(outxpos)
%         pxpwVec(1,:,itr) = px;
%         pxpwVec(2,:,itr) = w0p;
%     end
    
    % for each lighting direction, the viewing direction gets different
    % direction
    
    if(gpuEnable)
%         randPhase = exp(2*pi*1i*rand(1,1,1,1,1,pathsNum,'gpuArray'));
        randPhase = gpuArray(exp(2*pi*1i*rand(1,1,pathsNum)));
        
        us = feval(gpuFunc.singleScattering,us,lambda,gpuArray(x),(2*pi*randPhase./sqrt(px)), ...
            conv_vmf_l0.mu1, conv_vmf_l0.mu2, conv_vmf_l0.mu3, conv_vmf_l0.c, ...
            apertureVmf_v.mu1, apertureVmf_v.mu2, apertureVmf_v.mu3, apertureVmf_v.c,...
            apertureVmf_v.dir1,apertureVmf_v.dir2,apertureVmf_v.dir3);
    else
        randPhase = exp(2*pi*1i*rand(1,1,pathsNum));
        dz = cubeDist(x,box_min,box_max,dirv);
        throughputVmf_v = movmfThroughput(apertureVmf_v,lambda,x,1,sigt,dz);
        movmf_mult = movmfMultiple(conv_vmf_l0,throughputVmf_v,false);
        movmf_mult.alpha = movmf.alpha;
        [integrateMult] = movmfIntegrate(movmf_mult);
        us_res = integrateMult .* randPhase ./ sqrt(px);
        
        if(correlationFlag)
            us = us + squeeze(sum(us_res(1,1,:,:,:,:,:) .* conj(us_res(1,2,:,:,:,:,:)),1));
        else
            us = us + squeeze(sum(us_res(1,1,:,:,:,:),3));
        end
    end
   
    
    if(gpuEnable)
        activatedPaths = true(pathsNum,1,'gpuArray');
        el = feval(gpuFunc.integrateEl,el_zeros,w,w0p, ...
            throughputVmf_l.mu1, throughputVmf_l.mu2, throughputVmf_l.mu3, throughputVmf_l.c);
    else
        activatedPaths = true(pathsNum,1);
        rotatedMovmf_l = movmfRotate(movmf,w);
        throughputVmf_l_times_movmf = movmfMultiple(throughputVmf_l,rotatedMovmf_l,true);
        throughputVmf_l_times_movmf.c = throughputVmf_l_times_movmf.c - log(w0p);
        throughputVmf_l_times_movmf.alpha = movmf.alpha;
        el = movmfIntegrate(throughputVmf_l_times_movmf);
    end

    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
%     weight=albedo;
    

    
    % begin paths sampling loop
    while pL < maxBounces

        pL=pL+1;
        
%         if(outxpos)
%             xVec(:,pL,:,itr) = x;
%             wVec(:,pL,:,itr) = w;
%         end


        % calculate the complex volumetric throughput for the last
        % scattering event in case of multiple scattering
        if (pL>1)
                
            if(gpuEnable)
                randPhase = gpuArray(exp(2*pi*1i*rand(1,1,pathsNum)));
                
                u = feval(gpuFunc.multipleScattering,u, el, ...
                    lambda, x, ow, (2 * pi * randPhase ./ sqrt(px)), activatedPaths, ...
                    apertureVmf_v.mu1, apertureVmf_v.mu2, apertureVmf_v.mu3, apertureVmf_v.c, ...
                    apertureVmf_v.dir1,apertureVmf_v.dir2,apertureVmf_v.dir3);
            else
                randPhase = exp(2*pi*1i*rand(1,1,1,1,1,pathsNum));

                % rotate the scattering function on direction of ow
                rotatedMovmf_v = movmfRotate(movmf,ow);
                
                dz = cubeDist(x,box_min,box_max,dirv);
                throughputVmf_v = movmfThroughput(apertureVmf_v,lambda,x,1,sigt,dz);
                throughputVmf_v_times_movmf = movmfMultiple(throughputVmf_v,rotatedMovmf_v,true);
                throughputVmf_v_times_movmf.alpha = movmf.alpha;
                ev = movmfIntegrate(throughputVmf_v_times_movmf);
                u_res = (ev .* el ./ sqrt(px) .* randPhase);
                
                if(correlationFlag)
                    u = u + squeeze(sum(u_res(activatedPaths,1,:,:,:,:,:) .* conj(u_res(activatedPaths,2,:,:,:,:,:)),1));
                else
                    u = u + squeeze(sum(u_res(1,1,activatedPaths,:,:,:),3));
                end
            end
            
        end

       
        % advance to the next scattering event
        d=-log(-rand(1,1,pathsNum)+1)/(sigt);
        x = x+d.*w;

        % move to the next particle if the next scattering event is outside
        % the box
        outsideIdx = ...
            x(1,1,:) > box_max(1) | x(1,1,:) < box_min(1) | ...
            x(2,1,:) > box_max(2) | x(2,1,:) < box_min(2) | ...
            x(3,1,:) > box_max(3) | x(3,1,:) < box_min(3) ;
        
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
