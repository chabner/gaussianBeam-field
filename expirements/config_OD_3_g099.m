config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e3;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;
config.mcGpuV = true;
config.mcGpuL = true;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
% config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.scattgMFP = 100;
config.attMFP = 100;
config.boxDepth = 300;
config.boxAxial = 1e5;

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
config.focalPointsL.base = -20:1:20;
config.focalPointsL.xyGrid = true;
config.focalPointsL.plain = -config.boxDepth/2;
config.focalPointsL.dim = 3;

% -------------
% focal view points
config.focalPointsV.base = -40:1:40;
config.focalPointsV.xyGrid = true;
config.focalPointsV.plain = config.focalPointsL.plain;
config.focalPointsV.dim = 2;

% -------------
% focal illumination directions
config.focalDirectionsL.base = 0;
config.focalDirectionsL.xyGrid = false;
config.focalDirectionsL.dim = 4;

% -------------
% focal view directions
config.focalDirectionsV.base = 0;
config.focalDirectionsV.xyGrid = false;
config.focalDirectionsV.dim = 4;

% -------------
% shift v
config.focalPointsVshift.vShift = 0;
config.focalPointsVshift.dim = 4;

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.99; % HG parameter
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
% nf code, position:
% 1: random
% 2: exponent
% 3: gaussian multiple (not activated)
% 4: gaussian sum
%
% nf code, direction:
% 1: random
% 2: not activated
% 3: gaussian multiple (not activated)
% 4: gaussian sum
% 5: gaussian sum, same beam (only 54 is possible)

% ff code, position:
% 1: random
% 2: exponent

% ff code, direction:
% 1: random
% 2: sample from g0 (need to be provided)
% 3: sample from g0 - multimode

% choose 1 for random, and 3 for gaussian
config.sampleFlag = 54;

% kappa g parameter
config.kappaG = 100;
