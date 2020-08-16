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
config.MFP = 10;
config.boxDepth = 10;
config.boxAxial = 10;
config.boxShift = [0;0;0];

%% common parameters
config.parameters.lambda = [0.5,1,2];
config.wavelenght = @(lambda) lambda;

%% ff parameters
config.ff.parameters.v_x = -1:0.1:1;
config.ff.parameters.v_y = -1:0.1:1;
config.ff.parameters.l_x = -1:0.1:1;
config.ff.parameters.l_y = -1:0.1:1;

config.ff.v.x = @(v_x) v_x;
config.ff.v.y = @(v_y) v_y;
config.ff.v.z = @(v_x,v_y) sqrt(1 - v_x.^2 - v_y.^2);
config.ff.l.x = @(l_x) l_x;
config.ff.l.y = @(l_y) l_y;
config.ff.l.z = @(l_x,l_y) sqrt(1 - l_x.^2 - l_y.^2);

%% nf Parameters
config.nf.parameters.v_x = -5:0.125:5;
config.nf.parameters.v_y = -5:0.125:5;
config.nf.parameters.delta = 0:0.5:1.5;

%% Aperture

% mask of the gaussian lens
config.nf.mask_varL = 0.25;
config.nf.mask_varV = 0.25;

config.nf.focalPointsL.x = @(delta) -delta/2;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -5;

config.nf.focalPointsV.x = @(delta,v_x) -delta/2 + v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -5;

config.nf.focalDirectionsL.x = @() 0;
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @() 1;

config.nf.focalDirectionsV.x = @() 0;
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @() 1;

config.nf.focalPointsL.x_2 = @(delta) delta/2;
config.nf.focalPointsL.y_2 = @() 0;
config.nf.focalPointsL.z_2 = @() -5;

config.nf.focalPointsV.x_2 = @(delta,v_x) delta/2 + v_x;
config.nf.focalPointsV.y_2 = @(v_y) v_y;
config.nf.focalPointsV.z_2 = @() -5;

config.nf.focalDirectionsL.x_2 = @() 0;
config.nf.focalDirectionsL.y_2 = @() 0;
config.nf.focalDirectionsL.z_2 = @() 1;

config.nf.focalDirectionsV.x_2 = @() 0;
config.nf.focalDirectionsV.y_2 = @() 0;
config.nf.focalDirectionsV.z_2 = @() 1;

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
