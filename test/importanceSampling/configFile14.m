config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e0;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;
config.mcGpuV = true;
config.mcGpuL = false;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
% config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.scattgMFP = 10;
config.attMFP = 10;
config.boxDepth = 20;
config.boxAxial = 20;

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.05;
config.mask_varV = 0.05;

% -------------

% depth of lighting and viewa
config.focalPointsL_plain = -config.boxDepth/2 - 0;
config.focalPointsV_plain = -config.boxDepth/2 - 0;

% -------------
% focal illumination points
config.focalPointsL_base = 0:1:4;

% -------------
% focal view points
config.focalPointsV_base = -15:0.5:15;

% -------------
% focal directions
config.focalLDirections = 5 * deg2rad(0:1:4);
config.focalVDirections = 5 * deg2rad(0:1:4);

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
% the tens digit is for direction, and the first digit is for position
% choose 1 for random, and 3 for gaussian
config.sampleFlag = 14;

% kappa g parameter
config.kappaG = 100;