config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e3;
config.multiplePaths = 1e0;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.MFP = 10;
config.boxDepth = 20;
config.boxAxial = 25;
config.boxShift = [0;0;0];

%% ff parameters
config.ff.parameters.v_x = -0.7:0.04:0.7;
config.ff.parameters.v_y = -0.7:0.04:0.7;
config.ff.parameters.l_x = -0.7:0.04:0.7;
config.ff.parameters.l_y = -0.7:0.04:0.7;
config.ff.parameters.z_l_sign = [1,-1];
config.ff.parameters.z_v_sign = [1,-1];

config.ff.v.x = @(v_x) v_x;
config.ff.v.y = @(v_y) v_y;
config.ff.v.z = @(v_x,v_y,z_l_sign) z_l_sign .* sqrt(1 - v_x.^2 - v_y.^2);
config.ff.l.x = @(l_x) l_x;
config.ff.l.y = @(l_y) l_y;
config.ff.l.z = @(l_x,l_y,z_v_sign) z_v_sign .* sqrt(1 - l_x.^2 - l_y.^2);

%% nf Parameters
config.nf.parameters.v_x = -6:0.125:6;
config.nf.parameters.v_y = -6:0.125:6;
config.nf.parameters.delta = [-2:1:2,-2:1:2];

%% Aperture

% mask of the gaussian lens
config.nf.mask_varL = 0.25;
config.nf.mask_varV = 0.25;

config.nf.focalPointsL.x = @(delta) delta;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -10;

config.nf.focalPointsV.x = @(v_x) v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -10;

config.nf.focalDirectionsL.x = @(delta) reshape(sin(deg2rad([-2:1:2,-2:1:2])),size(delta));
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @(delta) reshape(cos(deg2rad([-2:1:2,-2:1:2])),size(delta));

config.nf.focalDirectionsV.x = @(delta) reshape(sin(deg2rad([-2:1:2,178:1:182])),size(delta));
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @(delta) reshape(cos(deg2rad([-2:1:2,178:1:182])),size(delta));

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.9; % HG parameter
config.forwardWeight = 1; % forward weight for hg scattering

% wmf mixture settings
% maximum number of mixtures
config.vmf_k = 11;
% maximal iterations
config.vmf_iterations = 1e5;
% samples in each axis
config.vmf_samples = 1e6;

%% importance sampling
config.sampleFlag = 11;
