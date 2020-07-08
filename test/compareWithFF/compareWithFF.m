%% backScattering_g02
clear
plotCompare('backScattering_g02.m',false);

%% backScattering_g09
plotCompare('backScattering_g09.m',false);

%% twoDirections
plotCompare('twoDirections.m',false);

%% largeAngles
plotCompare('largeAngles.m',true);

%% smallAngles
plotCompare('smallAngles.m',true);

function plotCompare(configName,onlyFull)
    clear config
    run(['expirements',filesep,configName]);
    [config] = preprocessConfig(config);
    
    tic                
    [u_nf,us_nf] = run_rendering(config);
    toc
    
    clear config;
    run(['expirements',filesep,configName]);
    
    % forward
    config.ff.parameters.z_l_sign = 1;
    config.ff.parameters.z_v_sign = 1;
    
    config = preprocessConfig_ff(config);

    tic
    [u_ff_forward,us_ff_forward] = run_farField(config);
    toc
        
    % backward
    config.ff.parameters.z_l_sign = -1;
    config.ff.parameters.z_v_sign = 1;
    
    config = preprocessConfig_ff(config);
    
    tic
    [u_ff_backward,us_ff_backward] = run_farField(config);
    toc
        
    % full
    config.ff.parameters.z_l_sign = [1,-1];
    config.ff.parameters.z_v_sign = [1,-1];
    
    config = preprocessConfig_ff(config);
    
    tic
    [u_ff_full,us_ff_full] = run_farField(config);
    toc
        
    Nl = size(u_nf,3);
    
    figure;

    if(~onlyFull)
        subplot(3,3,1)
        imagesc(abs([ ...
            reshape(u_nf, size(u_nf,1) , []) ; ...
            reshape(u_ff_full, size(u_ff_full,1) , []) ...
            ]));
       title('full')

        subplot(3,3,2)
        imagesc(abs([ ...
            reshape(u_nf(:,:,1:Nl/2), size(u_nf,1) , []) ; ...
            reshape(u_ff_full(:,:,1:Nl/2), size(u_ff_full,1) , []) ...
            ]));

        subplot(3,3,3)
        imagesc(abs([ ...
            reshape(u_nf(:,:,(Nl/2+1):end), size(u_nf,1) , []) ; ...
            reshape(u_ff_full(:,:,(Nl/2+1):end), size(u_ff_full,1) , []) ...
            ]));

        subplot(3,3,4)
        imagesc(abs([ ...
            reshape(u_nf, size(u_nf,1) , []) ; ...
            reshape(u_ff_forward, size(u_ff_forward,1) , []) ...
            ]));
       title('forward')

        subplot(3,3,5)
        imagesc(abs([ ...
            reshape(u_nf(:,:,1:Nl/2), size(u_nf,1) , []) ; ...
            reshape(u_ff_forward(:,:,1:Nl/2), size(u_ff_forward,1) , []) ...
            ]));

        subplot(3,3,6)
        imagesc(abs([ ...
            reshape(u_nf(:,:,(Nl/2+1):end), size(u_nf,1) , []) ; ...
            reshape(u_ff_forward(:,:,(Nl/2+1):end), size(u_ff_forward,1) , []) ...
            ]));

        subplot(3,3,7)
        imagesc(abs([ ...
            reshape(u_nf, size(u_nf,1) , []) ; ...
            reshape(u_ff_backward, size(u_ff_backward,1) , []) ...
            ]));
       title('backward')

       subplot(3,3,8)
        imagesc(abs([ ...
            reshape(u_nf(:,:,1:Nl/2), size(u_nf,1) , []) ; ...
            reshape(u_ff_backward(:,:,1:Nl/2), size(u_ff_backward,1) , []) ...
            ]));

        subplot(3,3,9)
        imagesc(abs([ ...
            reshape(u_nf(:,:,(Nl/2+1):end), size(u_nf,1) , []) ; ...
            reshape(u_ff_backward(:,:,(Nl/2+1):end), size(u_ff_backward,1) , []) ...
            ]));
    else
        imagesc(abs([ ...
            reshape(u_nf, size(u_nf,1) , []) ; ...
            reshape(u_ff_full, size(u_ff_full,1) , []) ...
            ]));
    end
end