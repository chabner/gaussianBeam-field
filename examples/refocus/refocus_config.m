%% Simulation
config.simulation.projectName = mfilename;
config.simulation.gpuNum = 0; % gpu number begins from 0
config.simulation.iterations = 10;  % Each iteration is x1024 paths
config.simulation.precision = 'double'; % single or double

%% Sample
config.medium.dim = 3;
config.medium.MFP = 2;
config.medium.boxDepth = 4;
config.medium.boxAxial = 10;

config.medium.boxShift = [0;0;0];

%% Refocus
config.refocus.sample_random = false; % false: tabulated
config.refocus.binary_aperture = false; % false: vMF (Gaussian) apodization

% at least one of the both must be true
config.refocus.sample_forward = true;
config.refocus.sample_backward = true;

% for maximal xy value 1 we get all directions, limited x and y values may 
% be considered for binary aperture
config.refocus.max_xy_value = 1;

% bias the attenuation to the central direction
config.refocus.bias_attenuation = true;

% in case of tabulated directions: should be 1/boxAxial to avoid aliasing
% we assume dl = dv = dldv
config.refocus.tabulated_dldv = 0.1;

% in case of random directions: random number of each aperture (for
% example, for 64 we get 64x4 directions)
config.refocus.random_directions_number = 128;

%% Aperture
config.aperture.mask_varL = 0.25;
config.aperture.mask_varV = 0.25;

% is the aperture normalized
config.aperture.is_normalized = true;

%% nf Parameters
config.nf.parameters.v_x = -4:0.5:4;
config.nf.parameters.v_y = -4:0.5:4;
config.nf.parameters.delta = [0,1];
config.nf.parameters.k = 2*pi;

config.nf.focalPointsL.x = @(delta) -delta/2;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -config.medium.boxDepth/2;

config.nf.focalPointsV.x = @(delta,v_x) -delta/2 + v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -config.medium.boxDepth/2;

config.nf.wavenumber.k = @(k) k;

config.nf.focalDirectionsL.x = @() 0;
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @() 1;

config.nf.focalDirectionsV.x = @() 0;
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @() 1;

config.nf.focalPointsL.x_2 = @(delta) delta/2;
config.nf.focalPointsL.y_2 = @() 0;
config.nf.focalPointsL.z_2 = @() -config.medium.boxDepth/2;

config.nf.focalPointsV.x_2 = @(delta,v_x) delta/2 + v_x;
config.nf.focalPointsV.y_2 = @(v_y) v_y;
config.nf.focalPointsV.z_2 = @() -config.medium.boxDepth/2;

config.nf.focalDirectionsL.x_2 = @() 0;
config.nf.focalDirectionsL.y_2 = @() 0;
config.nf.focalDirectionsL.z_2 = @() 1;

config.nf.focalDirectionsV.x_2 = @() 0;
config.nf.focalDirectionsV.y_2 = @() 0;
config.nf.focalDirectionsV.z_2 = @() 1;

config.nf.wavenumber.k_2 = @(k) k;

%% Scattering fnuction

% scattering type
% 1: random
% 2: from tabulated
% 3: HG
config.scatter.type = 3;
config.scatter.g = 0.2; % HG parameter

%% vmf mixture settings
% maximum number of mixtures
config.movmf.vmf_k = 5;
% maximal iterations
config.movmf.vmf_iterations = 1e5;
% samples in each axis
config.movmf.vmf_samples = 1e6;

%% importance sampling
% position sampling
% 1: random
% 2: exponential
% 3: unused
% 4: gaussian sum
config.sample.position_type = 1;

% direction sampling
% 1: random
% 2: unused
% 3: unused
% 4: gaussian sum
config.sample.direction_type = 1;

% sample the position and direction from the same beam
config.sample.same_beam = true;

% tabulation size (x1024 samples)
config.sample.z_sample_num = 100;
