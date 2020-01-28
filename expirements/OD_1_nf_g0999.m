config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e3;

% total correlations matrix
config.totalMatrix = 1000;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;
config.mcGpu = true;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.scattgMFP = 100;
config.attMFP = 100;
config.boxDepth = 100;
config.boxAxial = 1e4;

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------

% depth of lighting and viewa
config.focalPointsL_plain = -config.boxDepth/2 - 0;
config.focalPointsV_plain = -config.boxDepth/2 - 0;

% -------------
% focal illumination points
config.focalPointsL_base = [-2,2];

% -------------
% focal view points
config.focalPointsV_base = -6:0.125:6;

% -------------
% focal directions
config.focalLDirections = deg2rad(0:1:1);
config.focalVDirections = deg2rad(0:1:1);

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.999; % HG parameter
config.forwardWeight = 1; % forward weight for hg scattering

% wmf mixture settings
% maximum number of mixtures
config.vmf_k = 13;
% maximal iterations
config.vmf_iterations = 1e5;
% samples in each axis
config.vmf_samples = 1e6;

%% importance sampling
% the tens digit is for direction, and the first digit is for position
% choose 1 for random, and 3 for gaussian
config.sampleFlag = 44;

% kappa g parameter
config.kappaG = 100;