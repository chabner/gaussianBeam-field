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

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
config.focalPointsL.base = [-2:1:2,-2:1:2];
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
config.focalDirectionsL.theta = deg2rad([-2:1:2,-2:1:2]);
config.focalDirectionsL.phi = 0;
config.focalDirectionsL.dim = 2;

% -------------
% focal view directions
config.focalDirectionsV.theta = deg2rad([-2:1:2,178:1:182]);
config.focalDirectionsV.phi = 0;
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
config.sampleFlag = 11;
