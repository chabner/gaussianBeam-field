config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e2;
config.multiplePaths = 1e3;

% t/f if use gpu
config.useGpu = true;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
% config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.MFP = 10;
config.boxDepth = 20;
config.boxAxial = 20;

%% nf Parameters
config.nf.parameters.v_x = -15:0.5:15;
config.nf.parameters.v_y = -15:0.5:15;
config.nf.parameters.delta = 0:1:4;

%% Aperture

% mask of the gaussian lens
config.nf.mask_varL = 0.05;
config.nf.mask_varV = 0.05;

config.nf.focalPointsL.x = @(delta) -delta/2;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -10;

config.nf.focalPointsV.x = @(delta,v_x) -delta/2 + v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -10;

config.nf.focalDirectionsL.x = @(delta) reshape(sin(5*deg2rad(0:1:4)),size(delta));
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @(delta) reshape(cos(5*deg2rad(0:1:4)),size(delta));

config.nf.focalDirectionsV.x = @(delta) reshape(sin(5*deg2rad(0:1:4)),size(delta));
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @(delta) reshape(cos(5*deg2rad(0:1:4)),size(delta));

config.nf.focalPointsL.x_2 = @(delta) delta/2;
config.nf.focalPointsL.y_2 = @() 0;
config.nf.focalPointsL.z_2 = @() -10;

config.nf.focalPointsV.x_2 = @(delta,v_x) delta/2 + v_x;
config.nf.focalPointsV.y_2 = @(v_y) v_y;
config.nf.focalPointsV.z_2 = @() -10;

config.nf.focalDirectionsL.x_2 = @(delta) reshape(sin(5*deg2rad(0:1:4)),size(delta));
config.nf.focalDirectionsL.y_2 = @() 0;
config.nf.focalDirectionsL.z_2 = @(delta) reshape(cos(5*deg2rad(0:1:4)),size(delta));

config.nf.focalDirectionsV.x_2 = @(delta) reshape(sin(5*deg2rad(0:1:4)),size(delta));
config.nf.focalDirectionsV.y_2 = @() 0;
config.nf.focalDirectionsV.z_2 = @(delta) reshape(cos(5*deg2rad(0:1:4)),size(delta));

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.5; % HG parameter
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
