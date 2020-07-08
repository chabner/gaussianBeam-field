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

if(~isfield(config,'multiplePaths'))
    config.multiplePaths = 1;
end

ff_correlationActive = isfield(config.ff,'v.x_2');
refocus_active = isfield(config,'nf');

if(refocus_active)
    refocus_correlationActive = isfield(config.nf.focalPointsL,'x_2');
    config.nf.correlationActive = refocus_correlationActive;
else
    refocus_correlationActive = false;
end


if(ff_correlationActive && refocus_active)
    error('refocus only on field')
end

if(config.multiplePaths ~= 1 && refocus_active)
    error('refocus only with a single path')
end

if(~isfield(config.ff,'dir_v'))
    config.noAttDir = true;
end

%% Parameters
% ff parameters
paramsList = fieldnames(config.ff.parameters);
currDim = 4;

ff_dims = numel(paramsList);

for paramNum = 1:1:ff_dims
    param = config.ff.parameters.(paramsList{paramNum});
    param = param(:);
    permuteVec = 1:1:currDim;
    permuteVec(1) = currDim; permuteVec(end) = 1;
    param = permute(param,permuteVec);
    
    config.ff.parameters.(paramsList{paramNum}) = param;
    currDim = currDim + 1;
end

% refocus parameters
if(refocus_active)
    paramsList = fieldnames(config.nf.parameters);
    nf_dims = numel(paramsList);

    for paramNum = 1:1:numel(paramsList)
        param = config.nf.parameters.(paramsList{paramNum});
        param = param(:);
        permuteVec = 1:1:currDim;
        permuteVec(1) = currDim; permuteVec(end) = 1;
        param = permute(param,permuteVec);

        config.nf.parameters.(paramsList{paramNum}) = param;
        currDim = currDim + 1;
    end
    
    config.nf.ff_dims = ff_dims;
    config.nf.nf_dims = nf_dims;
end

%% Evaluate
config.v.x = eval(['config.ff.v.x',parseFunc(config.ff.v.x,'config.ff.parameters.')]);
config.v.y = eval(['config.ff.v.y',parseFunc(config.ff.v.y,'config.ff.parameters.')]);
config.v.z = eval(['config.ff.v.z',parseFunc(config.ff.v.z,'config.ff.parameters.')]);

config.l.x = eval(['config.ff.l.x',parseFunc(config.ff.l.x,'config.ff.parameters.')]);
config.l.y = eval(['config.ff.l.y',parseFunc(config.ff.l.y,'config.ff.parameters.')]);
config.l.z = eval(['config.ff.l.z',parseFunc(config.ff.l.z,'config.ff.parameters.')]);

if(~config.noAttDir)
    config.dir_v.x = eval(['config.ff.dir_v.x',parseFunc(config.ff.dir_v.x,'config.ff.parameters.')]);
    config.dir_v.y = eval(['config.ff.dir_v.y',parseFunc(config.ff.dir_v.y,'config.ff.parameters.')]);
    config.dir_v.z = eval(['config.ff.dir_v.z',parseFunc(config.ff.dir_v.z,'config.ff.parameters.')]);

    config.dir_l.x = eval(['config.ff.dir_l.x',parseFunc(config.ff.dir_l.x,'config.ff.parameters.')]);
    config.dir_l.y = eval(['config.ff.dir_l.y',parseFunc(config.ff.dir_l.y,'config.ff.parameters.')]);
    config.dir_l.z = eval(['config.ff.dir_l.z',parseFunc(config.ff.dir_l.z,'config.ff.parameters.')]);
end

if(ff_correlationActive)
    config.v.x_2 = eval(['config.ff.v.x_2',parseFunc(config.ff.v.x_2,'config.ff.parameters.')]);
    config.v.y_2 = eval(['config.ff.v.y_2',parseFunc(config.ff.v.y_2,'config.ff.parameters.')]);
    config.v.z_2 = eval(['config.ff.v.z_2',parseFunc(config.ff.v.z_2,'config.ff.parameters.')]);

    config.l.x_2 = eval(['config.ff.l.x_2',parseFunc(config.ff.l.x_2,'config.ff.parameters.')]);
    config.l.y_2 = eval(['config.ff.l.y_2',parseFunc(config.ff.l.y_2,'config.ff.parameters.')]);
    config.l.z_2 = eval(['config.ff.l.z_2',parseFunc(config.ff.l.z_2,'config.ff.parameters.')]);

    if(~config.noAttDir)
        config.dir_v.x_2 = eval(['config.ff.dir_v.x_2',parseFunc(config.ff.dir_v.x_2,'config.ff.parameters.')]);
        config.dir_v.y_2 = eval(['config.ff.dir_v.y_2',parseFunc(config.ff.dir_v.y_2,'config.ff.parameters.')]);
        config.dir_v.z_2 = eval(['config.ff.dir_v.z_2',parseFunc(config.ff.dir_v.z_2,'config.ff.parameters.')]);

        config.dir_l.x_2 = eval(['config.ff.dir_l.x_2',parseFunc(config.ff.dir_l.x_2,'config.ff.parameters.')]);
        config.dir_l.y_2 = eval(['config.ff.dir_l.y_2',parseFunc(config.ff.dir_l.y_2,'config.ff.parameters.')]);
        config.dir_l.z_2 = eval(['config.ff.dir_l.z_2',parseFunc(config.ff.dir_l.z_2,'config.ff.parameters.')]);
    end
    
    config.v.x = cat(2,config.v.x,config.v.x_2);
    config.v.y = cat(2,config.v.y,config.v.y_2);
    config.v.z = cat(2,config.v.z,config.v.z_2);
    
    config.l.x = cat(2,config.l.x,config.l.x_2);
    config.l.y = cat(2,config.l.y,config.l.y_2);
    config.l.z = cat(2,config.l.z,config.l.z_2);
    
    if(~config.noAttDir)
        config.dir_v.x = cat(2,config.dir_v.x,config.dir_v.x_2);
        config.dir_v.y = cat(2,config.dir_v.y,config.dir_v.y_2);
        config.dir_v.z = cat(2,config.dir_v.z,config.dir_v.z_2);

        config.dir_l.x = cat(2,config.dir_l.x,config.dir_l.x_2);
        config.dir_l.y = cat(2,config.dir_l.y,config.dir_l.y_2);
        config.dir_l.z = cat(2,config.dir_l.z,config.dir_l.z_2);
    end
end

if(refocus_active)
    config.nf.focalPointsL.eval.x = eval(['config.nf.focalPointsL.x',parseFunc(config.nf.focalPointsL.x)]);
    config.nf.focalPointsL.eval.y = eval(['config.nf.focalPointsL.y',parseFunc(config.nf.focalPointsL.y)]);
    config.nf.focalPointsL.eval.z = eval(['config.nf.focalPointsL.z',parseFunc(config.nf.focalPointsL.z)]);

    config.nf.focalPointsV.eval.x = eval(['config.nf.focalPointsV.x',parseFunc(config.nf.focalPointsV.x)]);
    config.nf.focalPointsV.eval.y = eval(['config.nf.focalPointsV.y',parseFunc(config.nf.focalPointsV.y)]);
    config.nf.focalPointsV.eval.z = eval(['config.nf.focalPointsV.z',parseFunc(config.nf.focalPointsV.z)]);

    config.nf.focalDirectionsL.eval.x = eval(['config.nf.focalDirectionsL.x',parseFunc(config.nf.focalDirectionsL.x)]);
    config.nf.focalDirectionsL.eval.y = eval(['config.nf.focalDirectionsL.y',parseFunc(config.nf.focalDirectionsL.y)]);
    config.nf.focalDirectionsL.eval.z = eval(['config.nf.focalDirectionsL.z',parseFunc(config.nf.focalDirectionsL.z)]);

    config.nf.focalDirectionsV.eval.x = eval(['config.nf.focalDirectionsV.x',parseFunc(config.nf.focalDirectionsV.x)]);
    config.nf.focalDirectionsV.eval.y = eval(['config.nf.focalDirectionsV.y',parseFunc(config.nf.focalDirectionsV.y)]);
    config.nf.focalDirectionsV.eval.z = eval(['config.nf.focalDirectionsV.z',parseFunc(config.nf.focalDirectionsV.z)]);
end

if(refocus_correlationActive)
    config.nf.focalPointsL.eval.x_2 = eval(['config.nf.focalPointsL.x_2',parseFunc(config.nf.focalPointsL.x_2)]);
    config.nf.focalPointsL.eval.y_2 = eval(['config.nf.focalPointsL.y_2',parseFunc(config.nf.focalPointsL.y_2)]);
    config.nf.focalPointsL.eval.z_2 = eval(['config.nf.focalPointsL.z_2',parseFunc(config.nf.focalPointsL.z_2)]);

    config.nf.focalPointsV.eval.x_2 = eval(['config.nf.focalPointsV.x_2',parseFunc(config.nf.focalPointsV.x_2)]);
    config.nf.focalPointsV.eval.y_2 = eval(['config.nf.focalPointsV.y_2',parseFunc(config.nf.focalPointsV.y_2)]);
    config.nf.focalPointsV.eval.z_2 = eval(['config.nf.focalPointsV.z_2',parseFunc(config.nf.focalPointsV.z_2)]);

    config.nf.focalDirectionsL.eval.x_2 = eval(['config.nf.focalDirectionsL.x_2',parseFunc(config.nf.focalDirectionsL.x_2)]);
    config.nf.focalDirectionsL.eval.y_2 = eval(['config.nf.focalDirectionsL.y_2',parseFunc(config.nf.focalDirectionsL.y_2)]);
    config.nf.focalDirectionsL.eval.z_2 = eval(['config.nf.focalDirectionsL.z_2',parseFunc(config.nf.focalDirectionsL.z_2)]);

    config.nf.focalDirectionsV.eval.x_2 = eval(['config.nf.focalDirectionsV.x_2',parseFunc(config.nf.focalDirectionsV.x_2)]);
    config.nf.focalDirectionsV.eval.y_2 = eval(['config.nf.focalDirectionsV.y_2',parseFunc(config.nf.focalDirectionsV.y_2)]);
    config.nf.focalDirectionsV.eval.z_2 = eval(['config.nf.focalDirectionsV.z_2',parseFunc(config.nf.focalDirectionsV.z_2)]);
    
    config.nf.focalPointsL.eval.x = cat(2,config.nf.focalPointsL.eval.x,config.nf.focalPointsL.eval.x_2);
    config.nf.focalPointsL.eval.y = cat(2,config.nf.focalPointsL.eval.y,config.nf.focalPointsL.eval.y_2);
    config.nf.focalPointsL.eval.z = cat(2,config.nf.focalPointsL.eval.z,config.nf.focalPointsL.eval.z_2);
    
    config.nf.focalPointsV.eval.x = cat(2,config.nf.focalPointsV.eval.x,config.nf.focalPointsV.eval.x_2);
    config.nf.focalPointsV.eval.y = cat(2,config.nf.focalPointsV.eval.y,config.nf.focalPointsV.eval.y_2);
    config.nf.focalPointsV.eval.z = cat(2,config.nf.focalPointsV.eval.z,config.nf.focalPointsV.eval.z_2);
    
    config.nf.focalDirectionsV.eval.x = cat(2,config.nf.focalDirectionsV.eval.x,config.nf.focalDirectionsV.eval.x_2);
    config.nf.focalDirectionsV.eval.y = cat(2,config.nf.focalDirectionsV.eval.y,config.nf.focalDirectionsV.eval.y_2);
    config.nf.focalDirectionsV.eval.z = cat(2,config.nf.focalDirectionsV.eval.z,config.nf.focalDirectionsV.eval.z_2);
    
    config.nf.focalDirectionsL.eval.x = cat(2,config.nf.focalDirectionsL.eval.x,config.nf.focalDirectionsL.eval.x_2);
    config.nf.focalDirectionsL.eval.y = cat(2,config.nf.focalDirectionsL.eval.y,config.nf.focalDirectionsL.eval.y_2);
    config.nf.focalDirectionsL.eval.z = cat(2,config.nf.focalDirectionsL.eval.z,config.nf.focalDirectionsL.eval.z_2);
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