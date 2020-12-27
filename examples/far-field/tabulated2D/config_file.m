config.simulation.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.simulation.gpuNum = 0; % gpu number begins from 0
config.simulation.iterations = 100;  % Each iteration is x1024 paths
config.simulation.precision = 'single'; % single or double
% config.simulation.cbs = ...; % is coherent back scattering is active

%% Sample
config.medium.dim = 2;
config.medium.MFP = 10;
config.medium.boxDepth = 30;
config.medium.boxAxial = 1e4;
config.medium.boxShift = [0;0];

%% ff parameters
config.ff.parameters.theta_v = 0:0.002:(2*pi);
config.ff.parameters.theta_l = [0,0.001,0.005,0.01];
config.ff.parameters.is_bw = [1,-1];

config.ff.v.x = @(theta_v,theta_l) sin(theta_v - theta_l/2);
config.ff.v.y = @(theta_v,theta_l) cos(theta_v - theta_l/2);

config.ff.l.x = @(theta_l,is_bw) sin(is_bw .* -theta_l/2);
config.ff.l.y = @(theta_l) cos(theta_l/2);

config.ff.wavenumber.k = @() 2*pi;

config.ff.v.x_2 = @(theta_v,theta_l) sin(theta_v + theta_l/2);
config.ff.v.y_2 = @(theta_v,theta_l) cos(theta_v + theta_l/2);

config.ff.l.x_2 = @(theta_l,is_bw) sin(is_bw .* theta_l/2);
config.ff.l.y_2 = @(theta_l) cos(theta_l/2);

config.ff.wavenumber.k_2 = @() 2*pi;

%% Scattering fnuction

% scattering type
% 1: random
% 2: tabulated
% 3: HG
config.scatter.type = 2;

% HG g parameter
% config.scatter.g = 0.5;

% Tabulated complex amplitude function
% config.scatter.f = ...

%% Sampling
% position sampling
% 1: random
% 2: exponential
config.sample.position_type = 2;

% direction sampling
% 1: random
% 2: from tabulated function
% 3: from hg function
config.sample.direction_type = 1;

% HG g parameter
% config.sample.g0 = 0.95;

% Tabulated complex amplitude function
% config.sample.f = tab_func.f;

% in case of correlation: if we sample from the average of two directions
config.sample.mean_l = true;
