%% Simulation
config.simulation.projectName = mfilename;
config.simulation.gpuNum = 0; % gpu number begins from 0
config.simulation.iterations = 100;  % Each iteration is x1024 paths
config.simulation.precision = 'single'; % single or double

%% Sample
config.medium.dim = 3;
config.medium.MFP = 10;
config.medium.boxDepth = 50;
config.medium.boxAxial = 1000;
config.medium.boxShift = [0;0;0];

%% Aperture
config.aperture.mask_varL = 0.25;
config.aperture.mask_varV = 0.25;

% is the aperture normalized
config.aperture.is_normalized = true;

%% nf Parameters
config.nf.parameters.v_x = -7:0.125:7;
config.nf.parameters.v_y = -7:0.125:7;
config.nf.parameters.delta = 0:0.25:1.5;
config.nf.parameters.k = 2*pi;

config.nf.focalPointsL.x = @(delta) delta;
config.nf.focalPointsL.y = @() 0;
config.nf.focalPointsL.z = @() -config.medium.boxDepth/2;

config.nf.focalPointsV.x = @(v_x) v_x;
config.nf.focalPointsV.y = @(v_y) v_y;
config.nf.focalPointsV.z = @() -config.medium.boxDepth/2;

config.nf.wavenumber.k = @(k) k;

config.nf.focalDirectionsL.x = @() 0;
config.nf.focalDirectionsL.y = @() 0;
config.nf.focalDirectionsL.z = @() 1;

config.nf.focalDirectionsV.x = @() 0;
config.nf.focalDirectionsV.y = @() 0;
config.nf.focalDirectionsV.z = @() 1;

%% Scattering fnuction

% scattering type
% 1: random
% 2: from tabulated
% 3: HG
config.scatter.type = 3;
config.scatter.g = 0.95; % HG parameter

%% vmf mixture settings
% maximum number of mixtures
config.movmf.vmf_k = 11;
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
config.sample.position_type = 4;

% direction sampling
% 1: random
% 2: unused
% 3: unused
% 4: gaussian sum
config.sample.direction_type = 4;

% sample the position and direction from the same beam
config.sample.same_beam = true;

% tabulation size (x1024 samples)
config.sample.z_sample_num = 10;

% min probability

% one of the two, or calculate the min probability for path
config.sample.test_rounds = 10; % each round is x1024 paths
config.sample.min_probability_percent = 0.01;

% or set the minimal probability
% config.sample.min_px0 = 0;
% config.sample.min_pw0 = 0;
