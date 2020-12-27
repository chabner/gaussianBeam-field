config.simulation.projectName = mfilename;

%% Rendering

% iterations to render single correlation image
config.simulation.gpuNum = 0; % gpu number begins from 0
config.simulation.iterations = 100;  % Each iteration is x1024 paths
config.simulation.precision = 'single'; % single or double
config.simulation.cbs = false; % is coherent back scattering is active

%% Sample
config.medium.dim = 3;
config.medium.MFP = 10;
config.medium.boxDepth = 30;
config.medium.boxAxial = 1e4;
config.medium.boxShift = [0;0;0];

%% ff parameters
config.ff.parameters.vx = -0.2:0.005:0.2;
config.ff.parameters.vy = -0.2:0.005:0.2;
config.ff.parameters.dl = [0,0.05,0.2,0.5];

config.ff.v.x = @(vx,dl) vx - dl/2;
config.ff.v.y = @(vy) vy;
config.ff.v.z = @(vx,vy,dl) sqrt(1 - (vx - dl/2).^2 - vy.^2);

config.ff.l.x = @(dl) -dl/2;
config.ff.l.y = @() 0;
config.ff.l.z = @(dl) sqrt(1 -(dl/2).^2);

config.ff.wavenumber.k = @() 2*pi;

config.ff.v.x_2 = @(vx,dl) vx + dl/2;
config.ff.v.y_2 = @(vy) vy;
config.ff.v.z_2 = @(vx,vy,dl) sqrt(1 - (vx + dl/2).^2 - vy.^2);

config.ff.l.x_2 = @(dl) dl/2;
config.ff.l.y_2 = @() 0;
config.ff.l.z_2 = @(dl) sqrt(1 -(dl/2).^2);

config.ff.wavenumber.k_2 = @() 2*pi;

%% Scattering fnuction

% scattering type
% 1: random
% 2: tabulated
% 3: HG
config.scatter.type = 3;

% HG g parameter
config.scatter.g = 0.9;

% Tabulated complex amplitude function
% config.scatter.f = 

%% Sampling
% position sampling
% 1: random
% 2: exponential
config.sample.position_type = 1;

% direction sampling
% 1: random
% 2: from tabulated function
% 3: from hg function
config.sample.direction_type = 1;

% HG g parameter
config.sample.g0 = 0.8;

% Tabulated complex amplitude function
% config.sample.f = 

% in case of correlation: if we sample from the average of two directions
config.sample.mean_l = true;
