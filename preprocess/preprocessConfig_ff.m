function [config] = preprocessConfig_ff(config)
%% Sample
config.box_min = [-config.boxAxial/2;-config.boxAxial/2;-config.boxDepth/2];
config.box_max = [config.boxAxial/2;config.boxAxial/2;config.boxDepth/2];

if(isfield(config,'boxShift'))
    config.box_min = config.box_min + config.boxShift;
    config.box_max = config.box_max + config.boxShift;
end

%% Scattering function
if(config.sctType == 3)
    config.sct_type = 3;
    config.ampfunc.g = config.g;
    config.ampfunc.forwardWeight = config.forwardWeight;
end

if(isfield(config,'g0'))
    config.sct_type = 3;
    config.ampfunc0.g = config.g0;
    config.ampfunc0.forwardWeight = config.forwardWeight;
else
    config.ampfunc0 = inf;
end
%% Importance sampling
% Sample 100 different samples, and take the median of px, in order to
% estimate bad samples
pxItersNum = 100;
minTorr = 1e-3;

box_w = config.box_max-config.box_min;
V=prod(box_w);

px = zeros(1,pxItersNum);

for iterNum = 1:1:pxItersNum
    switch mod(config.sampleFlag,10)
        case 1
            % uniform distribution
            px(iterNum) = 1;
        case 2
            % exponential distribution
            [~,tpx] = expSmpX(config.box_min,config.box_max,[0;0;1],1/config.attMFP);
            px(iterNum) = tpx;
    end
end

config.smpPreprocess(1).pxMin = median(px) * minTorr;

end