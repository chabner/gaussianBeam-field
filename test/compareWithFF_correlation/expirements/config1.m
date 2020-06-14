config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e3;
config.multiplePaths = 1e2;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
% config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.MFP = 10;
config.boxDepth = 10;
config.boxAxial = 10;
config.boxShift = [0;0;0];

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
delta = 0:0.5:1.5;
config.focalPointsL.dim = 3;

config.focalPointsL.xyGrid = false;
config.focalPointsL.base = -delta/2;
config.focalPointsL.plain = -config.boxDepth/2 - 0;

config.focalPointsL.base_2 = delta/2;
config.focalPointsL.plain_2 = -config.boxDepth/2 - 0;

% -------------
% focal view points
config.focalPointsV.xyGrid = true;
config.focalPointsV.dim = 2;

config.focalPointsV.base = -5:0.125:5;
config.focalPointsV.plain = -config.boxDepth/2 - 0;

config.focalPointsV.base_2 = -5:0.125:5;
config.focalPointsV.plain_2 = -config.boxDepth/2 - 0;

% -------------
% focal illumination directions
config.focalDirectionsL.xyGrid = false;
config.focalDirectionsL.dim = 3;

config.focalDirectionsL.base = 0 * config.focalPointsL.base;
config.focalDirectionsL.base_2 = 0 * config.focalPointsL.base;

% -------------
% focal view directions
config.focalDirectionsV.xyGrid = false;
config.focalDirectionsV.dim = 3;

config.focalDirectionsV.base = 0 * config.focalPointsL.base;
config.focalDirectionsV.base_2 = 0 * config.focalPointsL.base;

% -------------
% shift v
config.focalPointsVshift.vShift = config.focalPointsL.base;
config.focalPointsVshift.vShift_2 = config.focalPointsL.base_2;
config.focalPointsVshift.dim = 3;

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
