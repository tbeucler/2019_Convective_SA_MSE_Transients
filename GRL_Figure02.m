%% GRL_Figure02.m
% tbeucler - 6/17/2019
% Edited for Github repository of GRL submission
% tbeucler - 3/17/2019
% Spectra of MSE for all products

close all; fclose('all'); clearvars;

%% 1. Parameters

% Figure's appearance
fz = 16; % Fontsize
lw = 3; % Linewidth

% Subplot 2a parameters
% Simulations
LEG2a = {'0-5d CTRL','0-5d UNI-RAD','5-10d CTRL','5-10d UNI-RAD',...
    '40-80d CTRL','40-80d UNI-RAD'};
sim2a = {'LC_UNI_RAD','LC_CTRL'}; Nsim = length(sim2a);
TIT2a = ['$\mathrm{\left(a\right)\ \lambda^{-1}\times MSE\ Power',...
    '\ Spectrum\ for\ LC\ experiments}$'];
YLAB2a = '$\varphi_{H}/\lambda\ \left[\mathrm{kg^{2}\ m^{-4}}\right]$';
% Time periods to average MSE power spectrum for subplot 2a
COL2a = [[127 48 143];[218 89 33];[237 177 33]]/255; % Colors for different times
ls2a = {':','-'}; % Linestyles for different time periods
tmin_2a = [0 5 40];
tmax_2a = [5 10 80];
N2a = length(tmin_2a);

% Subplot 2b parameters
% Parameters that change for each dataset
COL2b = [[237 177 33];[237 177 33];[237 177 33];...
    [0 127.5 0];[0 127.5 0];[0 127.5 0];...
    [127.5 127.5 255];[0 0 0]]/255;
LEG2b = {'LC UNI-RAD','LC UNI-SEF','LC CTRL','NG UNI-RAD',...
    'NG UNI-SEF','NG CTRL','CERES','ERA'};
ls2b = {':','--','-',':','--','-','-','-'}; N2b = length(ls2b);
sim2b = {'LC_UNI_RAD_1day_av','LC_UNI_SEF_1day_av','LC_CTRL_1day_av',...
    'NG_UNI_RAD','NG_UNI_SEF','NG_CTRL','CERES','ERA5_1day_av'};
TIT2b = ['$\mathrm{\left(b\right)\ \lambda^{-1}\times MSE',...
    '\ Power\ Spectrum\ for\ all\ datasets}$'];
tmin_2b = [40 40 40 50 50 50 0 0];
tmax_2b = [80 80 80 80 80 80 1e4 1e4];
XLAB2b = '$\lambda\ \left[\mathrm{km}\right]$';
YLAB2b = '$\varphi_{H}/\lambda\ \left[\mathrm{kg^{2}\ m^{-4}}\right]$';
% Position adjustments
leg = 0.225; % Legend margin
marl = 0.1; % Left margin
marr = 2.5e-2; % Right margin
marv = 5e-2; % Vertical margin

% Physical constants
Lv = 2.5e6; % Latent heat of vaporization of water
spd = 24*3600; % Number of seconds per day

%% 2. Figure
figure
set(gcf,'Position',[50 50 1000 700]);

%% 2a. Subplot 2a
S(1) = subplot(2,1,1);
for iT = 1:N2a, tmin = tmin_2a(iT); tmax = tmax_2a(iT); % Time period
    for isim = 1:Nsim, load(['MAT_DATA',filesep,sim2a{isim}]);
        % Time coordinate
        t = DAT.t-DAT.t(1);
        [~,i1] = min(abs(t-tmin)); [~,i2] = min(abs(t-tmax)); TT = i1:i2;
        % Wavelength abscissa
        X = log10(DAT.lam_interp/1e3); Xmin = min(X(:));
        % Plot MSE spectrum normalized to variance units
        Y = log10(nanmean(DAT.Agg.mse(:,:,TT)/2.,3).*...
            nanmean(DAT.VAR.mse(:,TT)/Lv^2,2)./...
            (DAT.lam_interp'.*trapz(1./DAT.lam_interp,...
            nanmean(DAT.Agg.mse(:,:,TT)/2.,3))));
        P2a(isim,iT) = plot(X,Y,'Linewidth',lw,'color',COL2a(iT,:),...
            'Linestyle',ls2a{isim}); hold on; grid on;
        set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex');
    end
    % Subplot's label and appearance
    set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex',...
        'Xdir','reverse','XTickLabel','');
    
end; ylabel(YLAB2a,'Fontsize',fz,'Interpreter','Latex');
title(TIT2a,'Fontsize',fz,'Interpreter','Latex');
LEGa = legend([P2a(2,1) P2a(1,1) P2a(2,2) P2a(1,2) P2a(2,3) P2a(1,3)],...
    LEG2a,'Location','eastoutside','Interpreter','Latex','Fontsize',fz);

% Use logarithmic scale for y labels
G = gca; YTIK = G.YTickLabel;
for iy = 1:numel(YTIK), YTIK{iy}=strcat('$10^{',YTIK{iy},'}$'); end
set(gca,'YTickLabel',YTIK);

%% 2b. Subplot 2b
S(2) = subplot(2,1,2);
for isim = 1:N2b, load(['MAT_DATA',filesep,sim2b{isim}]);
    tmin = tmin_2b(isim); tmax = tmax_2b(isim); % Time period
    % Time coordinate
    t = DAT.t-DAT.t(1);
    [~,i1] = min(abs(t-tmin)); [~,i2] = min(abs(t-tmax)); TT = i1:i2;
    % Wavelength abscissa
    X = log10(DAT.lam_interp/1e3); Xmax = max(X(:));
    % Ordinate
    Y = log10(nanmean(DAT.Agg.mse(:,:,TT)/2,3).*...
        nanmean((Lv^(-2))*DAT.VAR.mse(:,TT),2)./...
        (DAT.lam_interp'.*trapz(1./...
        DAT.lam_interp,nanmean(DAT.Agg.mse(:,:,TT)/2,3)))); hold on;
    % Plot power MSE spectrum normalized to variance units
    P2b(isim) = plot(X,Y,'Linewidth',lw,'color',COL2b(isim,:),...
        'Linestyle',ls2b{isim}); hold on; grid on;
    set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex');
end; set(gca,'Xdir','reverse','Fontsize',fz); xlim([Xmin Xmax]);
% Axes and appearance of subplot 2b
ylabel(YLAB2b,'Fontsize',fz,'Interpreter','Latex');
xlabel(XLAB2b,'Fontsize',fz,'Interpreter','Latex');
title(TIT2b,'Fontsize',fz,'Interpreter','Latex');
LEGb = legend(P2b,LEG2b,'Location','Eastoutside',...
    'Interpreter','Latex','Fontsize',fz);

S(1).XLim = [Xmin Xmax]; % Set same x limits for subplot 2a

% Adjust position of subplots 2a,2b and corresponding legends
S(1).Position = [marl 0.5+marv 1-marl-marr-leg 0.5-2*marv];
LEGa.Position = [1-marr-leg 0.5+marv leg 0.5-2*marv];
S(2).Position = [marl 2*marv 1-marl-marr-leg 0.5-2*marv];
LEGb.Position = [1-marr-leg 2*marv leg 0.5-2*marv];

% Logarithmic scales for both axes of subplot 2b
G = gca; YTIK = G.YTickLabel; XTIK = G.XTickLabel;
for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
for iy = 1:numel(YTIK), YTIK{iy}=strcat('$10^{',YTIK{iy},'}$'); end
set(gca,'XTickLabel',XTIK,'YTickLabel',YTIK);

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'GRL_Fig02.pdf']); % Save plot