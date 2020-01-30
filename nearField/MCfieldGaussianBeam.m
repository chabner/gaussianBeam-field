function [u,us,um,u_nomean,xRep,wRep]=MCfieldGaussianBeam(sigt,albedo,box_min,box_max,xl,xv,signl,signv,varl,varv,dirl,dirv,maxItr,lambda,smpFlg,smpFunc,movmf,sct_type,ampfunc,ampfunc0)
%xl: a 3xNl or 2xNl vector of positions of where beam is focused (do not have to be
%on a gread)
%xv: 3xNv or 2xNv, center of viewing beam
%sign l, sign v : is beam facing upward (1) or downward (-1)
%varl, varv: scalars. beam std in fourier plane
u_nomean = 0;
xRep = 0;
wRep = 0;

% get the dimensions size
if((numel(box_max) ~= 2 && numel(box_max) ~= 3) || ...
        (size(box_max,2) ~= 1) || (any(size(box_max) ~= size(box_min))))
    error('Invalid box size');
end

dim = size(box_min,1);

unmeanl = [0;0;1];

Nl = size(xl,2);
Nv = size(xv,2);

%% Prepare for algorithm

% Initiate output parameters
u = zeros(Nv,Nl,class(xv));
us =zeros(Nv,Nl,class(xv));


% Box size
box_w = box_max-box_min;
V=prod(box_w);

apertureVmf_l = movmfAperture(varl,xl,-1,dirl);

apertureVmf_v = cell(1,Nl);
for lightNum = 1:1:Nl
    apertureVmf_v{lightNum} = movmfAperture(varv,xv,1,dirv(:,lightNum));
end

% threshold to begin kill particles with low weight
killThr=0.2;

%% Begin the main loop
for itr=1:maxItr
    if(mod(itr,1e3) == 0) 
        disp([num2str(round(1e2 * itr/maxItr)), '%'])
    end
    % Sample the first scatter
    % x: first scattering point
    % px: probability by which first point was sampled. Needed so that
    % importance sampling integration is weighted properly
    switch mod(smpFlg,10)
        case 1
            % uniform distribution
            x=rand(dim,1).*(box_w)+box_min;
            px=1;
        case 2
            % exponential distribution
            [x,px]=expSmpX(box_min,box_max,unmeanl,sigt(2));
        case 3
            [x,px,xMissIter] = smpVmfBeam(apertureVmf_l,smpFunc,box_min,box_max);
            px = px * V; % The reson I multiple with V - same prob as uniform
            xRep = xRep + xMissIter;
        case 4
            [x,px,xMissIter,n] = smpVmfBeamSum(apertureVmf_l,smpFunc,box_min,box_max);
            px = px * V; % The reson I multiple with V - same prob as uniform
            xRep = xRep + xMissIter;
    end
%      x
    px = max(px,smpFunc(1).pxMin);
    
    dz = cubeDist(x,box_min,box_max,-dirl);
    throughputVmf_l = movmfThroughput(apertureVmf_l,x,-signl,sigt(2),dz);
    conv_vmf_l0 = movmfConv(throughputVmf_l,movmf);
    
    % First scattering direction
    %implement an option to sample w isotropically (add apropriate input flag).
    %otherwise implement an option that samples one of the incoming
    %illumination of interest and sample a direction from that Gaussian. 
    %sampling probability w0p should be
    switch floor(smpFlg/10)
        case 1
            w=randn(dim,1); w=w/norm(w);
            w0p=1/sqrt(2^(dim-1)*pi);
        case 3
            [w,w0p] = vmfFirstDirection(conv_vmf_l0);
        case 4
            [w,w0p] = vmfFirstDirectionSum(conv_vmf_l0,n);
    end
    
    % for each lighting direction, the viewing direction gets different
    % direction
    
    randPhase = exp(2*pi*1i*rand);
    
    for lightNum = 1:1:Nl
        dz = cubeDist(x,box_min,box_max,dirv(:,lightNum));
        throughputVmf_v = movmfThroughput(apertureVmf_v{lightNum},x,signv,sigt(2),dz);
        movmf_mult = movmfMultiple(movmfPick(conv_vmf_l0,lightNum),throughputVmf_v,true,false);
        integrateMult = movmfIntegrate(movmf_mult);
        us(:,lightNum) = us(:,lightNum) + integrateMult * randPhase / sqrt(px);
    end
    
    rotatedMovmf_l = movmf;
    rotatedMovmf_l.mu = sign(rotatedMovmf_l.mu(:,:,3)) .* reshape(w,1,1,[]);
    
    throughputVmf_l_times_movmf = movmfMultiple(throughputVmf_l,rotatedMovmf_l,true);
    el = movmfIntegrate(throughputVmf_l_times_movmf) / w0p;

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
            
            randPhase = exp(2*pi*1i*rand);

            % rotate the scattering function on direction of ow
            rotatedMovmf_v = movmf;
            rotatedMovmf_v.mu = sign(rotatedMovmf_v.mu(:,:,3)) .* reshape(ow,1,1,[]);
            
            for lightNum = 1:1:Nl
                dz = cubeDist(x,box_min,box_max,dirv(:,lightNum));
                throughputVmf_v = movmfThroughput(apertureVmf_v{lightNum},x,signv,sigt(2),dz);
                throughputVmf_v_times_movmf = movmfMultiple(throughputVmf_v,rotatedMovmf_v,true);
                ev = movmfIntegrate(throughputVmf_v_times_movmf);
                u(:,lightNum) = u(:,lightNum)+ev(:)*el(lightNum)/sqrt(px) * randPhase;
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
u=u*sqrt(1/maxItr*V*sigt(2));
us=us*sqrt(1/maxItr*V*sigt(2));
um = u;

u = u+us;

end