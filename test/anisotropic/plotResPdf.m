%% Res tables
% line 1: vmf multiple scattering
%
% line 2: vmf single scattering
%
% line 3: ff multiple scattering
%
% line 4: ff single scattering
%
% column 1: intensity
%
% column 2: correlation - scaled to intensity
%
% column 3: correlation


clear

load('res.mat')

deltaIdx1 = 1;
deltaIdx2 = 2;



%% g - 098, OD 3 - dis 1
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 1
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 1
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;


deltaIdx1 = 1;
deltaIdx2 = 3;

%% g - 098, OD 3 - dis 2
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 2
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 2
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 4;

%% g - 098, OD 3 - dis 3
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 3
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 3
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 5;

%% g - 098, OD 3 - dis 4
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 4
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 4
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 6;

%% g - 098, OD 3 - dis 5
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 5
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 5
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 7;

%% g - 098, OD 3 - dis 6
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 6
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 6
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 8;

%% g - 098, OD 3 - dis 7
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 7
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 7
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 9;

%% g - 098, OD 3 - dis 8
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 8
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 8
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 10;

%% g - 098, OD 3 - dis 9
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 9
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 9
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;

deltaIdx1 = 1;
deltaIdx2 = 11;

%% g - 098, OD 3 - dis 10
s_vmf = vmf.config_OD_3_g098;
s_ff = ff.config_OD_3_g098;

plotRes;

%% g - 099, OD 3 - dis 10
s_vmf = vmf.config_OD_3_g099;
s_ff = ff.config_OD_3_g099;

plotRes;

%% g - 099, OD 6 - dis 10
s_vmf = vmf.config_OD_6_g099;
s_ff = ff.config_OD_6_g099;

plotRes;