function [u_small,u_big,us_small,us_big] = run_rendering(config) 

%% GB code
% config = buildConfig(config);

Nl = size(config.focalPointsL,2);
% Nv = numel(config.focalPointsV_base);

if(config.mcGpu)    
%     for lightNum = 1:1:Nl
%         config.focalPointsL{lightNum} = gpuArray(config.focalPointsL{lightNum});
        config.focalPointsV = gpuArray(config.focalPointsV);
%         config.focalDirectionsL{lightNum} = gpuArray(config.focalDirectionsL{lightNum});
        config.focalDirectionsV = gpuArray(config.focalDirectionsV);
%         config.box_min = gpuArray(config.box_min);
%         config.box_max = gpuArray(config.box_max);
%     end
end

for corrIter = 1:1:1
%     parfor lightNum = 1:1:Nl
%     for lightNum = 1:1:Nl
        [u,us] = ...
            MCfieldGaussianBeam( ...
          [1/config.scattgMFP,1/config.attMFP] ,           ... sigt
          1,                       ... albedo
          config.box_min,                 ... box_min
          config.box_max,                 ... box_max
          config.focalPointsL,            ... xl
          config.focalPointsV,            ... xv
          1,                       ... signl
          1,                       ... signv
          config.mask_varL,        ... varl
          config.mask_varV,        ... varv
          config.focalDirectionsL, ... dirl
          config.focalDirectionsV, ... dirv
          config.iterationsRender, ... maxItr
          config.wavelenght,       ... lambda
          config.sampleFlag,                       ... smpFlg
          config.smpPreprocess,        ... smpFunc
          config.movmf,                       ... movmf
          config.sct_type,                       ... sct_type
          config.ampfunc,                 ... ampfunc
          config.ampfunc0           ... ampfunc0
       );
%     end
end

if(config.mcGpu)
    u = gather(u);
    us = gather(us);
end

u_small = u(1:(numel(config.smallGrid)^2),:);
u_small = reshape(u_small,numel(config.smallGrid),numel(config.smallGrid),Nl);

u_big = u(((numel(config.smallGrid)^2)+1):end,:);
u_big = reshape(u_big,numel(config.bigGrid),numel(config.bigGrid),Nl);

us_small = us(1:(numel(config.smallGrid)^2),:);
us_small = reshape(us_small,numel(config.smallGrid),numel(config.smallGrid),Nl);

us_big = us(((numel(config.smallGrid)^2)+1):end,:);
us_big = reshape(us_big,numel(config.bigGrid),numel(config.bigGrid),Nl);

end