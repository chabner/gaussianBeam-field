clear

run(['test',filesep,'compareWithFF',filesep,'configFile.m']);
config = preprocessConfig(config);
[u_nf,us_nf] = run_rendering(config);

l = -0.5:0.02:0.5;
v = -0.5:0.02:0.5;
[u_ff,us_ff] = farField_render(config,l,v);


figure, imagesc(abs( reshape(u_nf, size(u_nf,1) , []) ));
figure, imagesc(abs( reshape(u_ff, size(u_ff,1) , []) ));

figure, imagesc(abs([ ...
    reshape(u_nf, size(u_nf,1) , []) ; ...
    reshape(u_ff, size(u_ff,1) , []) ...
    ]));