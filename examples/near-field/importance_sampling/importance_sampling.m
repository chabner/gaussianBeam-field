clear
addpath(genpath(['..',filesep,'..',filesep,'..',filesep]));
paths_num = [1,10,100,1000];
Np = numel(paths_num);

figure

%% random only

for pNum = 1:1:Np

    subplot(3,Np,pNum)
    run('config_file.m');
    config.sample.position_type = 1;
    config.sample.direction_type = 1;
    config.simulation.iterations = paths_num(pNum);
    config = preprocess_near_field(config);
    tic
    imagesc(config.nf.parameters.v_x,config.nf.parameters.v_y,abs(run_near_field(config)));
    toc
    
    colorbar
    title([num2str(paths_num(pNum) * 1024),' itetaions'])
    clear config

    if(pNum == 1)
        ylabel('random sampling')
    end
end

%% position only
for pNum = 1:1:numel(paths_num)

    subplot(3,Np,pNum + Np)
    run('config_file.m');
    config.sample.position_type = 4;
    config.sample.direction_type = 1;
    config.simulation.iterations = paths_num(pNum);
    config = preprocess_near_field(config);
    tic
    imagesc(config.nf.parameters.v_x,config.nf.parameters.v_y,abs(run_near_field(config)));
    toc
    
    colorbar
    title([num2str(paths_num(pNum) * 1024),' itetaions'])
    clear config

    if(pNum == 1)
        ylabel('position sampling')
    end
end

%% direction + position
for pNum = 1:1:numel(paths_num)

    subplot(3,Np,pNum + 2*Np)
    run('config_file.m');
    config.sample.position_type = 4;
    config.sample.direction_type = 4;
    config.simulation.iterations = paths_num(pNum);
    config = preprocess_near_field(config);
    tic
    imagesc(config.nf.parameters.v_x,config.nf.parameters.v_y,abs(run_near_field(config)));
    toc
    
    colorbar
    title([num2str(paths_num(pNum) * 1024),' itetaions'])
    clear config

    if(pNum == 1)
        ylabel('position and direction sampling')
    end
end
