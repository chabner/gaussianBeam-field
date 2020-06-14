function [u,us,u_nomean,xRep,xVec,wVec,pxpwVec]=MCfieldMOvMF2D(sigt,albedo,box_min,box_max,maxItr,lambda,movmf,sct_type,ampfunc,smpFlag,smpFunc,...
    apertureVmf_l_1,dirl_1,apertureVmf_v_1,dirv_1, ...
    apertureVmf_l_2,dirl_2,apertureVmf_v_2,dirv_2)
u_nomean = 0;
xRep = 0;
movmf.dim = movmf.dim(1:4);

outxpos = (nargout > 4);
maxBounces = 100;

if(outxpos)
    xVec = inf(2,maxBounces,1,maxItr);
    wVec = inf(2,maxBounces,1,maxItr);
    pxpwVec = inf(2,1,maxItr);
end

correlationFlag = (nargin > 15);

% get the dimensions size
if((numel(box_max) ~= 2 && numel(box_max) ~= 3) || ...
        (size(box_max,2) ~= 1) || (any(size(box_max) ~= size(box_min))))
    error('Invalid box size');
end

dim = size(box_min,1);

%% Prepare for algorithm

% Initiate output parameters
u_size = max(apertureVmf_l_1.dim,apertureVmf_v_1.dim);
u_size = u_size(2:end);

u = complex(zeros(u_size));
us = u;

% Box size
box_w = box_max-box_min;
V=prod(box_w);

% threshold to begin kill particles with low weight
killThr=0.2;

if(correlationFlag)
    apertureVmf_l = movmfUnite2D(apertureVmf_l_1, apertureVmf_l_2);
else
    apertureVmf_l = apertureVmf_l_1;
end


%% Begin the main loop
for itr=1:maxItr
    if(mod(itr,1e3) == 0) 
        disp([num2str(round(1e2 * itr/maxItr)), '%'])
    end
    switch mod(smpFlag,10)
        case 1
            % Sample the first scatter
            % x: first scattering point
            % uniform distribution
            x = rand(dim,1).*(box_w)+box_min;
            px = 1;
        case 4
            [x,px,xMissIter,n] = smpVmfBeamSum2D(apertureVmf_l,smpFunc,box_min,box_max);
            px = px * V; % The reson I multiple with V - same prob as uniform
            xRep = xRep + xMissIter;
    end
%     px = max(px,smpFunc(1).pxMin);
    dz = cubeDist(x,box_min,box_max,-dirl_1);
    throughputVmf_l_1 = movmfThroughput2D(apertureVmf_l_1,x,-1,sigt(2),dz);
%     conv_vmf_l0_1_1 = movmfConv2D2(throughputVmf_l_1,movmf,dirv_1,apertureVmf_v_1.dirDim);
    conv_vmf_l0_1 = movmfConv2D3(throughputVmf_l_1,movmf,dirv_1,apertureVmf_v_1.dirDim);
    
    if(correlationFlag)
        dz = cubeDist(x,box_min,box_max,-dirl_2);
        throughputVmf_l_2 = movmfThroughput2D(apertureVmf_l_2,x,-1,sigt(2),dz);
        conv_vmf_l0_2 = movmfConv2D3(throughputVmf_l_2,movmf,dirv_2,apertureVmf_v_2.dirDim);
    end
    
    % First scattering direction
    %implement an option to sample w isotropically (add apropriate input flag).
    %otherwise implement an option that samples one of the incoming
    %illumination of interest and sample a direction from that Gaussian. 
    %sampling probability w0p should be
    switch floor(smpFlag/10)
        case 1
            w=randn(dim,1); w=w/norm(w);
            w0p=1/sqrt(2^(dim-1)*pi);
        case 5
            if(correlationFlag)
                conv_vmf_l0 = movmfUnite2D(conv_vmf_l0_1, conv_vmf_l0_2);
                [w,w0p] = vmfFirstDirectionSum_sameBeam2D(conv_vmf_l0,n);
            else
                [w,w0p] = vmfFirstDirectionSum_sameBeam2D(conv_vmf_l0_1,n);
            end         
    end
    
    if(outxpos)
        pxpwVec(1,1,itr) = px;
        pxpwVec(2,1,itr) = w0p;
    end
    
    % for each lighting direction, the viewing direction gets different
    % direction
    randPhase = exp(2*pi*1i*rand);

    dz = cubeDist(x,box_min,box_max,dirv_1);
    throughputVmf_v_1 = movmfThroughput2D(apertureVmf_v_1,x,1,sigt(2),dz);
    movmf_mult = movmfMultiple2D(conv_vmf_l0_1,throughputVmf_v_1,false);
    integrateMult = movmfIntegrate2D(movmf_mult);
    us_1 = squeeze(integrateMult * randPhase / sqrt(px));
    
    if(correlationFlag)
        dz = cubeDist(x,box_min,box_max,dirv_2);
        throughputVmf_v_2 = movmfThroughput2D(apertureVmf_v_2,x,1,sigt(2),dz);
        movmf_mult = movmfMultiple2D(conv_vmf_l0_2,throughputVmf_v_2,false);
        integrateMult = movmfIntegrate2D(movmf_mult);
        us_2 = squeeze(integrateMult * randPhase / sqrt(px));
        
        u = u + us_1 .* conj(us_2);
        us = us + us_1 .* conj(us_2);
    else
        u = u + us_1;
        us = us + us_1;
    end
    
    rotatedMovmf_l = movmfRotate2D(movmf,w);
    
    throughputVmf_l_times_movmf = movmfMultiple2D(throughputVmf_l_1,rotatedMovmf_l,true);
    throughputVmf_l_times_movmf.c = throughputVmf_l_times_movmf.c - log(w0p);
    el_1 = movmfIntegrate2D(throughputVmf_l_times_movmf);
    
    if(correlationFlag)
        throughputVmf_l_times_movmf = movmfMultiple2D(throughputVmf_l_2,rotatedMovmf_l,true);
        throughputVmf_l_times_movmf.c = throughputVmf_l_times_movmf.c - log(w0p);
        el_2 = movmfIntegrate2D(throughputVmf_l_times_movmf);
    end

    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
    weight=albedo;
    
    % begin paths sampling loop
    while pL < maxBounces

        pL=pL+1;
        
        if(outxpos)
            xVec(:,pL,1,itr) = x;
            wVec(:,pL,1,itr) = w;
        end

        % calculate the complex volumetric throughput for the last
        % scattering event in case of multiple scattering
        if (pL>1)
            
            randPhase = exp(2*pi*1i*rand);

            % rotate the scattering function on direction of ow
            rotatedMovmf_v = movmfRotate2D(movmf,ow);
                
            dz = cubeDist(x,box_min,box_max,dirv_1);
            throughputVmf_v_1 = movmfThroughput2D(apertureVmf_v_1,x,1,sigt(2),dz);
            throughputVmf_v_times_movmf = movmfMultiple2D(throughputVmf_v_1,rotatedMovmf_v,true);
            ev = movmfIntegrate2D(throughputVmf_v_times_movmf);
            u_1 = squeeze(ev .* el_1 * randPhase / sqrt(px));
            
            if(correlationFlag)
                dz = cubeDist(x,box_min,box_max,dirv_2);
                throughputVmf_v_2 = movmfThroughput2D(apertureVmf_v_2,x,1,sigt(2),dz);
                throughputVmf_v_times_movmf = movmfMultiple2D(throughputVmf_v_2,rotatedMovmf_v,true);
                ev = movmfIntegrate2D(throughputVmf_v_times_movmf);
                u_2 = squeeze(ev .* el_2 * randPhase / sqrt(px));
                
                u = u + u_1 .* conj(u_2);
            else
                u = u + u_1;
            end
        end

       
        % advance to the next scattering event
        d=-log(-rand+1)/(sigt(1));
        x = x+d*w;

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
        w=smpampfunc_general(ow, sct_type,ampfunc);
    end
  
end

%% Normalization
if(correlationFlag)
    u=u*(1/(maxItr + xRep)*V*sigt(2));
    us=us*(1/(maxItr + xRep)*V*sigt(2));    
else
    u=u*sqrt(1/(maxItr + xRep)*V*sigt(2));
    us=us*sqrt(1/(maxItr + xRep)*V*sigt(2));
end



end
