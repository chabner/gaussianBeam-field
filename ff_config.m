config.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.iterationsRender = 4e3;
config.multiplePaths = 1e0;

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
config.boxDepth = 10;
config.boxAxial = 10;
config.boxShift = [0;0;0];

%% ff Parameters

config.ff.parameters.v_x = -1:0.1:1;
config.ff.parameters.v_y = -1:0.1:1;
config.ff.parameters.l_x = -1:0.1:1;
config.ff.parameters.l_y = -1:0.1:1;

%% ff Directions

% ff directions
config.ff.v.x = @(v_x) v_x;
config.ff.v.y = @(v_y) v_y;
config.ff.v.z = @(v_x,v_y) sqrt(1 - v_x.^2 - v_y.^2);

config.ff.l.x = @(l_x) l_x;
config.ff.l.y = @(l_y) l_y;
config.ff.l.z = @(l_x,l_y) sqrt(1 - l_x.^2 - l_y.^2);

% ff directions for attenuation
config.ff.dir_v.x = @(v_x) v_x;
config.ff.dir_v.y = @(v_y) v_y;
config.ff.dir_v.z = @(v_x,v_y) sqrt(1 - v_x.^2 - v_y.^2);

config.ff.dir_l.x = @(l_x) l_x;
config.ff.dir_l.y = @(l_y) l_y;
config.ff.dir_l.z = @(l_x,l_y) sqrt(1 - l_x.^2 - l_y.^2);

%% nf Parameters
config.nf.parameters.v_x = -5:0.2:5;
config.nf.parameters.v_y = -5:0.2:5;
config.nf.parameters.delta = 0:1:10;

%% Aperture

% mask of the gaussian lens
config.nf.mask_varL = 0.25;
config.nf.mask_varV = 0.25;

config.nf.focalPointsL.x = @(delta) -delta/2;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -500;

config.nf.focalPointsL.x_2 = @(delta) delta/2;
config.nf.focalPointsL.y_2 = @() 0;
config.nf.focalPointsL.z_2 = @() -500;

config.nf.focalPointsV.x = @(delta,v_x) -delta/2 + v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -500;

config.nf.focalPointsV.x_2 = @(delta,v_x) delta/2 + v_x;
config.nf.focalPointsV.y_2 = @(v_y) v_y;
config.nf.focalPointsV.z_2 = @() -500;

config.nf.focalDirectionsL.x = @(delta) (delta/2) / 250;
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @(delta) sqrt(1 - ((delta/2) / 250).^2);

config.nf.focalDirectionsL.x_2 = @(delta) (-delta/2) / 250;
config.nf.focalDirectionsL.y_2 = @() 0;
config.nf.focalDirectionsL.z_2 = @(delta) sqrt(1 - ((delta/2) / 250).^2);

config.nf.focalDirectionsV.x = @(delta) (delta/2) / 250;
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @(delta) sqrt(1 - ((delta/2) / 250).^2);

config.nf.focalDirectionsV.x_2 = @(delta) (-delta/2) / 250;
config.nf.focalDirectionsV.y_2 = @() 0;
config.nf.focalDirectionsV.z_2 = @(delta) sqrt(1 - ((delta/2) / 250).^2);

%% Scattering fnuction

% scattering type (only hg is implemented right now)
% HG
config.sctType = 3;
config.g = 0.97; % HG parameter
config.forwardWeight = 1; % forward weight for hg scattering

%% importance sampling
config.sampleFlag = 11;
