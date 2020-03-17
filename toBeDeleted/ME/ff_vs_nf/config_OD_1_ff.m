% vAngles = deg2rad(-10:0.1:10)/5;

farFactor = 1e0;

config.boxDepth = 100;
config.focalPointsV.base = (-3000:50:3000)*farFactor;
lBase = ([0:40:200,300:100:1000,2000:1000:4000])*farFactor;
config.focalPointsL.base = [-lBase/2,lBase/2];
config.focalPointsL.plain = -config.boxDepth/2 - 20000*farFactor;

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
config.scattgMFP = 100;
config.attMFP = 100;
config.boxAxial = 1e2;

%% Aperture

% mask of the gaussian lens
config.mask_varL = 1e2;
config.mask_varV = 1e2;

% -------------
% focal illumination points
config.focalPointsL.xyGrid = false;
config.focalPointsL.dim = 2;

% -------------
% focal view points
config.focalPointsV.xyGrid = true;
config.focalPointsV.plain = config.focalPointsL.plain;
config.focalPointsV.dim = 3;

% -------------
% focal illumination directions
config.focalDirectionsL.base = 0 * config.focalPointsL.base;
config.focalDirectionsL.xyGrid = false;
config.focalDirectionsL.dim = 2;

% -------------
% focal view directions
config.focalDirectionsV.base = 0 * config.focalPointsL.base;
config.focalDirectionsV.xyGrid = false;
config.focalDirectionsV.dim = 2;

% -------------
% shift v
config.focalPointsVshift.vShift = config.focalPointsL.base;
config.focalPointsVshift.dim = 2;

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.97; % HG parameter
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
config.sampleFlag = 32;

% kappa g parameter
config.kappaG = 100;

% g0 for ff
config.g0 = 0.97;