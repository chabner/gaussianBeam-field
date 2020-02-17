config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e2;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;
config.mcGpuV = true;
config.mcGpuL = true;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
% config.rng = 1;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.scattgMFP = 100;
config.attMFP = 100;
config.boxDepth = 100;
config.boxAxial = 1e4;
config.boxShift = [1.2;0.1;2.3];

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
config.focalPointsL.base = -50:1:50;
config.focalPointsL.xyGrid = true;
config.focalPointsL.plain = -50:10:50;
config.focalPointsL.dim = 2;

% -------------
% focal view points
config.focalPointsV.base = -50:1:50;
config.focalPointsV.xyGrid = true;
config.focalPointsV.plain = -50:10:50;
config.focalPointsV.dim = 2;

% -------------
% focal illumination directions
config.focalDirectionsL.base = deg2rad(-1:1:1);
config.focalDirectionsL.xyGrid = true;
config.focalDirectionsL.dim = 3;

% -------------
% focal view directions
config.focalDirectionsV.base = deg2rad(-1:1:1);
config.focalDirectionsV.xyGrid = true;
config.focalDirectionsV.dim = 4;

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
config.sampleFlag = 54;

% kappa g parameter
config.kappaG = 100;