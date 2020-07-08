config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 1e0;

% t/f if use gpu (for fitting algorithm)
config.useGpu = true;

% set rng before rendering.
% comment to avoid seeting rng.
% possible to rng 'shuffle'
% config.rng = 5;

%% Sample
config.dimNum = 3;
config.wavelenght = 1;
config.MFP = 5e2;
config.boxDepth = 1e4;
config.boxAxial = 1e8;

%% nf Parameters
config.nf.parameters.v_x = [0,100];
config.nf.parameters.v_z = [0,100];
config.nf.parameters.theta = [deg2rad([-10,0,10]),deg2rad(180 + [-10,0,10])];

%% Aperture

% mask of the gaussian lens
config.nf.mask_varL = 0.25;
config.nf.mask_varV = 0.25;

config.nf.focalPointsL.x = @() 0;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() 0;

config.nf.focalPointsV.x = @(v_x) v_x;
config.nf.focalPointsV.y = @() 0;
config.nf.focalPointsV.z = @(v_z) v_z;

config.nf.focalDirectionsL.x = @(theta) sin(theta);
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @(theta) cos(theta);

config.nf.focalDirectionsV.x = @(theta) sin(theta);
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @(theta) cos(theta);

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
config.sampleFlag = 54;
