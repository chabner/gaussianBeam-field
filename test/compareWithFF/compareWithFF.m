clear

run(['test',filesep,'compareWithFF',filesep,'configFile.m']);
config = preprocessConfig(config);
[u_nf,us_nf] = run_rendering(config);

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