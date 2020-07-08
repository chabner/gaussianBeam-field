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


%% nf Parameters
config.nf.parameters.v_x = (-40:1:40)/4;
config.nf.parameters.v_y = (-40:1:40)/4;
config.nf.parameters.delta = 0:1:3;

%% Aperture

% mask of the gaussian lens
config.nf.mask_varL = 0.25;
config.nf.mask_varV = 0.25;

config.nf.focalPointsL.x = @(delta) -delta/2;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -150;

config.nf.focalPointsV.x = @(delta,v_x) -delta/2 + v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -150;

config.nf.focalDirectionsL.x = @() 0;
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @() 1;

config.nf.focalDirectionsV.x = @() 0;
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @() 1;

config.nf.focalPointsL.x_2 = @(delta) delta/2;
config.nf.focalPointsL.y_2 = @() 0;
config.nf.focalPointsL.z_2 = @() -150;

config.nf.focalPointsV.x_2 = @(delta,v_x) delta/2 + v_x;
config.nf.focalPointsV.y_2 = @(v_y) v_y;
config.nf.focalPointsV.z_2 = @() -150;

config.nf.focalDirectionsL.x_2 = @() 0;
config.nf.focalDirectionsL.y_2 = @() 0;
config.nf.focalDirectionsL.z_2 = @() 1;

config.nf.focalDirectionsV.x_2 = @() 0;
config.nf.focalDirectionsV.y_2 = @() 0;
config.nf.focalDirectionsV.z_2 = @() 1;

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

