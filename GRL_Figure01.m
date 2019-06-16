%% GRL_Figure01.m
% tbeucler - 6/15/2019
% Figure01 of Github repository of GRL submission
% tbeucler - 3/13/2019
% Snapshots of moist static energy anomaly

close all; fclose('all'); clearvars;

%% 1. Parameters and initialization

% Figure parameters
fz = 10; % Fontsize
LIM_kgm2 = [-12 12]; % Colorbar limits in kg/m2
LIM_MJm2 = [2635 2825]; % Colorbar limits in MJ/m2
lw = 5; % Linewidth
marg = 5e-2; % Left margin
Xfig = 1000;  % Width of figure
Yfig = 700;  % Height of figure
asp = Xfig/Yfig; % Aspect ratio
Xsub = (1-2*marg); % Subplot's abscissa
Ysub = Xsub*asp/18; % Subplot's ordinate

% Physical constants
erad = 6371e3; % Earth radius
deg = pi/180; % Degree in radian
Lv = 2.5e6; % Latent heat of vaporization of water

% Name of dataset and subplot's position
TIT_array = {'(a) ERA: Total MSE on Jan 1, 2010',...
    '(b) ERA: Time-averaged MSE from Jan 1, 2010 to Dec 31, 2014',...
    '(c) ERA: Transient MSE on Jan 1, 2010',...
    '(d) NG CTRL: Transient MSE for the control experiment',...
    '(e) NG UNI-RAD: Transient MSE with horizontally-uniform radiation',...
    '(f) NG UNI-SEF: Transient MSE with horizontally-uniform surface fluxes',...
    '(g) LC CTRL','(h) LC UNI-RAD','(i) LC UNI-SEF'};
NAM_array = {'ERA_Jan_1','ERA_Time_mean','ERA_Transient','NG_CTRL',...
    'NG_UNI_RAD','NG_UNI_SEF','LC_CTRL','LC_UNI_RAD','LC_UNI_SEF'};
MIN_array = [1 4 7 10 13 16 19 20 21];
MAX_array = [3 6 9 12 15 18 19 20 21];

% Load data
load('MAT_DATA/Fig01.mat');

%% 2. Plot moist static energy fields as color maps
figure
set(gcf,'Position',[50 50 Xfig Yfig]);

for isub = 1:9, S(isub) = subplot(7,3,MIN_array(isub):MAX_array(isub));
    NAM = NAM_array{isub}; % Name of dataset
    % Define what to plot
    if isub<=2, TOPLOT = Fig01.(NAM).MSE'/1e6; % Convert to MJ/m2
        COL = BLUE; LIM = LIM_MJm2;
    else, TOPLOT = (Fig01.(NAM).MSE'-mean(mean(Fig01.(NAM).MSE)))/Lv; % Convert to kg/m2
        COL = REDBLUE; LIM = LIM_kgm2;
    end
    % Plot it
    P(isub) = pcolor(Fig01.(NAM).x,Fig01.(NAM).y,TOPLOT);
    set(P(isub),'LineStyle','none'); colormap(gca,COL); % With no black lines
    C(isub) = colorbar; caxis(LIM); % Colorbar
    title(TIT_array(isub),'Fontsize',fz,'Interpreter','Latex');
    set(gca,'TickLabelInterpreter','Latex','Fontsize',fz,'YTick','');
    if isub~=9||isub~=6, set(gca,'XTickLabel',''); end
    if isub>1&&isub<9, set(C(isub),'Visible','off'); end
end

% Set top colorbar
set(C(1),'Orientation','horizontal','Position',[marg 1-2*marg 1-2*marg marg/2],...
    'TickLabelInterpreter','Latex','Fontsize',fz);
set(C(1).Label,'String','$\mathrm{Total\ MSE\ \left[MJ\ m^{-2}\right]}$',...
    'Interpreter','Latex','Fontsize',1.5*fz);
% Set bottom colorbar
set(C(9),'Orientation','horizontal','Position',[marg 1.5*marg 1-2*marg marg/2],...
    'TickLabelInterpreter','Latex','Fontsize',fz);
set(C(9).Label,'String','$\mathrm{Transient\ MSE/L_{v}\ \left[kg\ m^{-2}\right]}$',...
    'Interpreter','Latex','Fontsize',1.5*fz,'Position',[0 -1.5 0]);

% Adjust subplot's positions
for isub = 1:6, S(isub).Position = [marg Xsub-1.5*isub*Ysub Xsub Ysub];
end
for isub = 7:9
    S(isub).Position = [marg+Xsub*(isub-7)/3 ...
        Xsub-10.5*Ysub+Ysub*(12.4-5)/(2*12.4) Xsub/3 Ysub*5/12.4];
    set(S(isub),'Linewidth',lw);
end

%  Add scale at the bottom of the long channel's plot
A1 = annotation('doublearrow','Linewidth',lw/2,...
    'X',[S(8).Position(1) S(8).Position(1)+S(8).Position(3)],...
    'Y',[S(8).Position(2)-marg/3 S(8).Position(2)-marg/3],'Units','normalized');
A2 = annotation('textbox','String','12,288 km','Interpreter','Latex','Linestyle','none',...
    'Fontsize',1.33*fz,'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Position',[S(8).Position(1) S(8).Position(2)-marg S(8).Position(3) marg/2]);

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'GRL_Fig01.pdf']); % Save plot