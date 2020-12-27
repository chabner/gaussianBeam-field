%% Simulation
config.simulation.projectName = mfilename;
config.simulation.gpuNum = 0; % gpu number begins from 0
config.simulation.iterations = 20;  % Each iteration is x1024 paths
config.simulation.precision = 'single'; % single or double

%% Sample
config.medium.dim = 3;
config.medium.MFP = 2;
config.medium.boxDepth = 4;
config.medium.boxAxial = 10;
config.medium.boxShift = [0;0;0];

%% Aperture
config.aperture.mask_varL = 0.25;
config.aperture.mask_varV = 0.25;

% is the aperture normalized
config.aperture.is_normalized = true;

%% nf Parameters
config.nf.parameters.k = 2*pi;

config.nf.focalPointsL.x = @() 0;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() 0;

config.nf.focalPointsV.x = @() 0;
config.nf.focalPointsV.y = @() 0;
config.nf.focalPointsV.z = @() 0;

config.nf.wavenumber.k = @(k) k;

config.nf.focalDirectionsL.x = @() 0;
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @() 1;

config.nf.focalDirectionsV.x = @() 0;
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @() 1;

config.nf.focalPointsL.x_2 = @() 0;
config.nf.focalPointsL.y_2 = @() 0;
config.nf.focalPointsL.z_2 = @() 0;

config.nf.focalPointsV.x_2 = @() 0;
config.nf.focalPointsV.y_2 = @() 0;
config.nf.focalPointsV.z_2 = @() 0;

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
config.scatter.g = 0.24; % HG parameter

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
