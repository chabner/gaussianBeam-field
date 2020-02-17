clear

% run(['test',filesep,'compareWithFF',filesep,'configFile.m']);
run('configFile.m');
[config,gpuFunc] = preprocessConfig(config);

Nv = numel(config.focalPointsV.base);
Nl = numel(config.focalPointsL.base);
                
[u_nf,us_nf] = run_rendering(config,gpuFunc);

u_nf = reshape(u_nf,Nl,Nv,Nv);
us_nf = reshape(us_nf,Nl,Nv,Nv);
u_nf = permute(u_nf,[2,3,1]);
us_nf = permute(us_nf,[2,3,1]);

l = -1:0.02:1;
v = -1:0.02:1;
[u_ff,us_ff] = farField_render(config,l,v);


figure, imagesc(abs( reshape(u_nf, size(u_nf,1) , []) ));
figure, imagesc(abs( reshape(u_ff, size(u_ff,1) , []) ));

figure, imagesc(abs([ ...
    reshape(u_nf, size(u_nf,1) , []) ; ...
    reshape(u_ff, size(u_ff,1) , []) ...
    ]));

figure, imagesc(abs([ ...
    reshape(us_nf, size(u_nf,1) , []) ; ...
    reshape(us_ff, size(u_ff,1) , []) ...
    ]));