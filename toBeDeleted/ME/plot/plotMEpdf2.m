% close all
clear

load('ff')

OD1color = [0, 0.4470, 0.7410];
OD5color = [0.8500, 0.3250, 0.0980];
OD10color = [0.9290, 0.6940, 0.1250];

f = figure;
f.Position = [0,0,450,400];

hold on
plotDisplacement(expirementsRes.config_OD_1_ff.config, expirementsRes.config_OD_1_ff.C, OD1color,':');
plotDisplacement(expirementsRes.config_OD_10_ff.config, expirementsRes.config_OD_10_ff.C, OD10color,':');

xlabel('\Delta [\mum]')
saveFigure(f,'me_displacement');

%% ME - angle

f = figure;
f.Position = [0,0,450,400];
hold on

config = expirementsRes.config_OD_1_ff.config;
DL = abs(config.focalPointsL.plain) + abs(config.boxDepth/2);
angles =  atan(config.focalPointsL.base / DL);
plotAngle(expirementsRes.config_OD_1_ff.C,angles,OD1color,':');

load('MEcor_Jan192020_od1.mat')
angles = atan(l(1,:));
plotAngle2(ncC, angles,'m',':');

config = expirementsRes.config_OD_10_ff.config;
DL = abs(config.focalPointsL.plain) + abs(config.boxDepth/2);
angles =  atan(config.focalPointsL.base / DL);
plotAngle(expirementsRes.config_OD_10_ff.C,angles,OD10color,':');

load('MEcor_Jan192020_od10.mat')
angles = atan(l(1,:));
plotAngle2(ncC, angles,'g',':');

a1 = gca;
legend_f = figure;
legend_f.Position = [0,0,450,400];
copyobj(a1,legend_f);

figure(f)
xlabel('\theta [deg]')
saveFigure(f,'me_angle');

% figure(legend_f)
% saveLegend(legend_f,'me_legend');



function plotDisplacement(config, C, specifyColor, specifyLine, ~)
    currXaxis = config.focalPointsL.base;

    C = abs(C).^2;
    C(C > 1e4) = 0;
    C(isnan(C)) = 0;
    A = squeeze(sum(sum(C)));
    
    if(nargin == 6)
        A = A([1,3:end]);
        currXaxis = currXaxis([1,3:end]);
    end
    
    plot(currXaxis/2, A/A(1),'Color',specifyColor,'LineStyle',specifyLine,'lineWidth',2);
    
    xlim([0,40])
    ylim([0,1])
end

function plotAngle(C, angles, specifyColor, specifyLine, ~)
    
    C = abs(C).^2;
    C(C > 1e4) = 0;
    C(isnan(C)) = 0;
    A = squeeze(max(max(C)));
    
    if(nargin == 6)
        A = A([1,3:end]);
        angles = angles([1,3:end]);
    end
    
    plot(rad2deg(angles), A/A(1),'Color',specifyColor,'LineStyle',specifyLine,'lineWidth',2);
    
    xlim([0,10])
    ylim([0,1.01])
end

function plotAngle2(C, angles, specifyColor, specifyLine)
    
    plot(rad2deg(angles), C,'Color',specifyColor,'LineStyle',specifyLine,'lineWidth',2);
    
    xlim([0,10])
    ylim([0,1.01])
end

function saveFigure(f,saveName)
    set(gca, 'FontName', 'Times New Roman')
    set(gca,'FontSize',20);

%     yticks([])

    % make no edges
    ax = gca;
    box on
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    
    ax_width = ax_width - 0.01;
    left = left + 0.005;

    ax_height = ax_height - 0.01;
    bottom = bottom + 0.005;
    
    ax.Position = [left bottom ax_width ax_height];
    
    save2pdf(['pdf',filesep,saveName,'.pdf'],f,600);
    
    close(f);
end

function saveLegend(f,saveName)
    set(gca, 'FontName', 'Times New Roman')
    set(gca,'FontSize',14);

    yticks([])
    xticks([])

    % make no edges
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    
    ax_width = ax_width - 0.01;
    left = left + 0.005;

    ax_height = ax_height - 0.01;
    bottom =2;
    
    lgd = columnlegend(3,{ ...
        ['OD 01',newline,'z^o = 0'],['OD 01',newline,'z^o = L'],['OD 01',newline,'Far Field'], ...
        ['OD 05',newline,'z^o = 0'],['OD 05',newline,'z^o = L'],['OD 05',newline,'Far Field'], ...
        ['OD 10',newline,'z^o = 0'],['OD 10',newline,'z^o = L'],['OD 10',newline,'Far Field']}, ...
        'Location','southoutside', 'padding', 0.3);
    
    lgd.Position(1) = 0.05;
    lgd.Position(2) = -0.35;
    lgd.Position(3) = 0.9;
    lgd.Position(4) = 0;
    
    ax.Position = [left bottom ax_width ax_height];
    
    save2pdf(['pdf',filesep,saveName,'.pdf'],f,600,2.3);
    
    close(f);
end