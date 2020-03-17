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
config.scattgMFP = 1e8;
config.attMFP = 1e8;
config.boxDepth = 1e4;
config.boxAxial = 1e8;

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
config.focalPointsL.base = 0;
config.focalPointsL.xyGrid = false;
config.focalPointsL.plain = 0;
config.focalPointsL.dim = 2;

% -------------
% focal view points
config.focalPointsV.base = [0,100];
config.focalPointsV.xyGrid = false;
config.focalPointsV.plain = [0,100];
config.focalPointsV.dim = 2;

% -------------
% focal illumination directions
config.focalDirectionsL.base = deg2rad([-80,0,80]);
config.focalDirectionsL.xyGrid = false;
config.focalDirectionsL.dim = 3;

% -------------
% focal view directions
config.focalDirectionsV.base = deg2rad([-80,0,80]);
config.focalDirectionsV.xyGrid = false;
config.focalDirectionsV.dim = 3;


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
config.sampleFlag = 44;

% kappa g parameter
config.kappaG = 100;