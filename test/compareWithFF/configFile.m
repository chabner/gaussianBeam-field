config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e3;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;
config.mcGpuV = false;
config.mcGpuL = false;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.scattgMFP = 10;
config.attMFP = 10;
config.boxDepth = 20;
config.boxAxial = 20;
config.boxShift = [1.2;0.1;2.3];

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
config.focalPointsL.base = -5:1:5;
config.focalPointsL.xyGrid = false;
config.focalPointsL.plain = -config.boxDepth/2 - 0;
config.focalPointsL.dim = 2;

% -------------
% focal view points
config.focalPointsV.base = -6:0.125:6;
config.focalPointsV.xyGrid = true;
config.focalPointsV.plain = -config.boxDepth/2 - 0;
config.focalPointsV.dim = 3;

% -------------
% focal illumination directions
config.focalDirectionsL.base = 5 * deg2rad(-5:1:5);
config.focalDirectionsL.xyGrid = false;
config.focalDirectionsL.dim = 2;

% -------------
% focal view directions
config.focalDirectionsV.base = 5 * deg2rad(5:-1:-5);
config.focalDirectionsV.xyGrid = false;
config.focalDirectionsV.dim = 2;

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.2; % HG parameter
config.forwardWeight = 1; % forward weight for hg scattering

% wmf mixture settings
% maximum number of mixtures
config.vmf_k = 5;
% maximal iterations
config.vmf_iterations = 1e5;
% samples in each axis
config.vmf_samples = 1e6;

%% importance sampling
% the tens digit is for direction, and the first digit is for position
% choose 1 for random, and 3 for gaussian
config.sampleFlag = 11;

% kappa g parameter
config.kappaG = 100;