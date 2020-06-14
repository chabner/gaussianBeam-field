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
config.MFP = 100;
config.boxDepth = 300;
config.boxAxial = 3000;

%% Aperture

% mask of the gaussian lens
config.mask_varL = 0.25;
config.mask_varV = 0.25;

% -------------
% focal illumination points
config.focalPointsL.xyGrid = false;
config.focalPointsL.dim = 3;

config.focalPointsL.base = -(0:1:3)/2;
config.focalPointsL.plain = -config.boxDepth/2;

config.focalPointsL.base_2 = (0:1:3)/2;
config.focalPointsL.plain_2 = -config.boxDepth/2;

% -------------
% focal view points
config.focalPointsV.xyGrid = true;
config.focalPointsV.dim = 2;

config.focalPointsV.base = (-40:1:40)/4;
config.focalPointsV.plain = config.focalPointsL.plain;

config.focalPointsV.base_2 = (-40:1:40)/4;
config.focalPointsV.plain_2 = config.focalPointsL.plain_2;

% -------------
% focal illumination directions
config.focalDirectionsL.theta = 0;
config.focalDirectionsL.phi = 0;
config.focalDirectionsL.dim = 4;

config.focalDirectionsL.theta_2 = 0;
config.focalDirectionsL.phi_2 = 0;

% -------------
% focal view directions
config.focalDirectionsV.theta = 0;
config.focalDirectionsV.phi = 0;
config.focalDirectionsV.dim = 4;

config.focalDirectionsV.theta_2 = 0;
config.focalDirectionsV.phi_2 = 0;

% -------------
% shift v
config.focalPointsVshift.vShift = config.focalPointsL.base;
config.focalPointsVshift.vShift_2 = config.focalPointsL.base_2;

config.focalPointsVshift.dim = 3;

%% Scattering fnuction

% scattering type
% 2: tabulated. angel is 0 to pi in 3d, and 0 to 2pi in 2D.
% config.sctType = 2;
% config.pdf = [0,1,2,3]; % pdf is amplitude function. That is, sqrt of phase function.

% 3: HG
% config.sctType = 3;
% config.g = 0.99; % HG parameter
% config.forwardWeight = 1; % forward weight for hg scattering

% wmf mixture settings
% maximum number of mixtures
config.vmf_k = 7;

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

