function [config] = preprocessConfig_ff(config)
%% Sample
if(config.dimNum == 3)
    config.box_min = [-config.boxAxial/2;-config.boxAxial/2;-config.boxDepth/2];
    config.box_max = [config.boxAxial/2;config.boxAxial/2;config.boxDepth/2];
end

if(config.dimNum == 2)
    config.box_min = [-config.boxAxial/2;-config.boxDepth/2];
    config.box_max = [config.boxAxial/2;config.boxDepth/2];
end

if(isfield(config,'boxShift'))
    config.box_min = config.box_min + config.boxShift;
    config.box_max = config.box_max + config.boxShift;
end

if(isfield(config,'multiplePaths'))
    config.iterationsRender = config.iterationsRender * config.multiplePaths;
end

%% Scattering function
if(config.sctType == 3)
    config.sct_type = 3;
    config.ampfunc.g = config.g;
    config.ampfunc.forwardWeight = config.forwardWeight;
end

if(config.sctType == 2)
    config.sct_type = 2;
    
    if(config.dimNum == 3)
        theta = linspace(0,pi,numel(config.pdf));
        f = config.pdf;
        
        config.ampfunc.cs = mean((abs(f).^2).*sin(theta))*2*pi^2;
        config.ampfunc.samplePdf = (abs(f(:).').^2) .* sin(theta) ./ ...
            sum((abs(f(:).').^2) .* sin(theta));
        config.ampfunc.sampleCdf = cumsum(config.ampfunc.samplePdf);

        config.ampfunc.evalPdf = (abs(f(:).').^2) ./ config.ampfunc.cs;
        config.ampfunc.evalAmp = config.ampfunc.evalPdf .^ 0.5;
        config.ampfunc.sampleIcdf = invCDFvec(config.ampfunc.sampleCdf);
    end
    
    if(config.dimNum == 2)
        f = config.pdf;
        
        config.ampfunc.cs = mean(abs(f).^2)*2*pi;
        config.ampfunc.samplePdf = (abs(f(:).').^2) ./ sum((abs(f(:).').^2));
        config.ampfunc.sampleCdf = cumsum(config.ampfunc.samplePdf);

        config.ampfunc.evalPdf = (abs(f(:).').^2) ./ config.ampfunc.cs;
%         config.ampfunc.evalAmp = (config.ampfunc.evalPdf .^ 0.5) .* exp(1i*angle(f(:).'));
        config.ampfunc.evalAmp = complex(config.ampfunc.evalPdf .^ 0.5);
        config.ampfunc.sampleIcdf = invCDFvec(config.ampfunc.sampleCdf);
    end
end

if(config.sampleFlag ~= 66)
    if(isfield(config,'g0'))
        config.sct_type = 3;
        config.ampfunc0.g = config.g0;
        config.ampfunc0.forwardWeight = config.forwardWeight;
    else
        config.ampfunc0 = inf;
    end
end

%% Importance sampling
% Sample 100 different samples, and take the median of px, in order to
% estimate bad samples
pxItersNum = 100;
minTorr = 1e-3;

box_w = config.box_max-config.box_min;
V=prod(box_w);

px = zeros(1,pxItersNum);
d0 = zeros(config.dimNum,1);
d0(end) = 1;

for iterNum = 1:1:pxItersNum
    switch mod(config.sampleFlag,10)
        case 1
            % uniform distribution
            px(iterNum) = 1;
        case 2
            % exponential distribution
            [~,tpx] = expSmpX(config.box_min,config.box_max,d0,1/config.MFP);
            px(iterNum) = tpx;
        case 6
            px(iterNum) = inf;
    end
end

config.smpPreprocess(1).pxMin = median(px) * minTorr;

end