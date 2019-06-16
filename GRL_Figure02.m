%% GRL_Figure02.m
% tbeucler - 3/17/2019
% Spectra of MSE for all products

close all; fclose('all'); clearvars;

fz = 16; lw = 3; % Fontsize and Linewidth
COL = [[127 48 143];[218 89 33];[237 177 33]]/255; % MSE power spectra at different times
cmap = [0 0 0; 117 112 179; 217 150 0; 231 41 138; 27 158 119]/255; % LW,SW,SEF,ADV colors

av6 = 24;
spd = 24*3600;
Lv = 2.5e6;

figure
set(gcf,'Position',[50 50 1000 700]);

%% Subplot (a)
S(1) = subplot(2,1,1);
for iT = 1:3
    
    if iT==1, t1 = 0; t2 = 5;
    elseif iT==2, t1 = 5; t2 = 10;
    elseif iT==3, t1 = 50; t2 = 80;
    end
    
    for isim = [1 3]
        if isim<3, load(['MAT_DATA/UCP2_2212109_SST',num2str(isim),'_rad','cam'],'DAT');
        elseif isim==3, load(['MAT_DATA/UCP2_2212109_SST',num2str(300),'_rad','cam'],'DAT');
        end
        if isim==1, ls = ':';
        elseif isim==2, ls = '--';
        elseif isim==3, ls = '-';
        end
        t = DAT.t-DAT.t(1);
        [~,i1] = min(abs(t-t1)); [~,i2] = min(abs(t-t2)); TT = i1:i2;
        
        X = log10(DAT.lam_interp/1e3); Xmin = min(X(:));
        Y = log10(nanmean(DAT.Agg.mse(:,:,TT)/2.,3).*...
            nanmean(DAT.VAR(:,TT)/Lv^2,2)./...
            (DAT.lam_interp'.*trapz(1./DAT.lam_interp,nanmean(DAT.Agg.mse(:,:,TT)/2.,3))));
        P0(isim,iT) = plot(X,Y,'Linewidth',lw,'color',COL(iT,:),'Linestyle',ls); hold on; grid on;
        set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex');
    end
    
    ylabel('$\varphi_{H}/\lambda\ \left[\mathrm{kg^{2}\ m^{-4}}\right]$','Fontsize',fz,'Interpreter','Latex');
    set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex',...
        'Xdir','reverse','XTickLabel','');
    
end; title('$\mathrm{\left(a\right)\ \lambda^{-1}\times MSE\ Power\ Spectrum\ for\ LC\ experiments}$',...
    'Fontsize',fz,'Interpreter','Latex');
LEG1 = legend([P0(3,1) P0(1,1) P0(3,2) P0(1,2) P0(3,3) P0(1,3)],...
    {'0-5d CTRL','0-5d UNI-RAD','5-10d CTRL','5-10d UNI-RAD','50-80d CTRL','50-80d UNI-RAD'},...
    'Location','eastoutside','Interpreter','Latex','Fontsize',fz);

G = gca; YTIK = G.YTickLabel;
for iy = 1:numel(YTIK), YTIK{iy}=strcat('$10^{',YTIK{iy},'}$'); end
set(gca,'YTickLabel',YTIK);

%% Subplot (b)
S(2) = subplot(2,1,2);
% Plot Long-Channel
t1 = 50; t2 = 80;
for isim = 1:3
    if av6>0, strav6 = ['_av6',num2str(av6)]; else, strav6 = ''; end
    if isim<3, load(['MAT_DATA/UCP2_2212109_SST',num2str(isim),'_rad','cam',strav6],'DAT');
    elseif isim==3, load(['MAT_DATA/UCP2_2212109_SST',num2str(300),'_rad','cam',strav6],'DAT');
    end
    if isim==1, ls = ':';
    elseif isim==2, ls = '--';
    elseif isim==3, ls = '-';
    end
    t = DAT.t-DAT.t(1);
    [~,i1] = min(abs(t-t1)); [~,i2] = min(abs(t-t2)); TT = i1:i2;
    
    X = log10(DAT.lam_interp/1e3);
    Y = log10(nanmean(DAT.Agg.mse(:,:,TT)/2.,3).*...
        nanmean(DAT.VAR(:,TT)/Lv^2,2)./...
        (DAT.lam_interp'.*trapz(1./DAT.lam_interp,nanmean(DAT.Agg.mse(:,:,TT)/2.,3))));
    P1(isim) = plot(X,Y,'Linewidth',lw,'color',COL(iT,:),'Linestyle',ls); hold on; grid on;
    set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex');
end
% Plot NG
for isim = 1:3
    if isim<3, load(['MAT_DATA/GRL_GetNG_',num2str(isim+1)],'DAT');
    elseif isim==3, load('MAT_DATA/GRL_GetNG_1','DAT');
    end
    if isim==1, ls = ':';
    elseif isim==2, ls = '--';
    elseif isim==3, ls = '-';
    end
    t = DAT.t-DAT.t(1);
    [~,i1] = min(abs(t-t1)); [~,i2] = min(abs(t-t2)); TT = i1:i2;
    
    X = log10(DAT.lam_interp/1e3);
    Y = log10(nanmean(DAT.Agg.mse(:,:,TT)/2.,3).*...
        nanmean(DAT.VAR.mse(:,TT)/Lv^2,2)./...
        (DAT.lam_interp.*trapz(1./DAT.lam_interp,nanmean(DAT.Agg.mse(:,:,TT)/2.,3))));
    P2(isim) = plot(X,Y,'Linewidth',lw,'color',[0 0.5 0],'Linestyle',ls); hold on; grid on;
    set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex','Xdir','reverse');
end



% Plot CERES
load('C:\Users\Tom\Desktop\Radiative_convective_instability\Cloud_rad_water_vap\MAT_DATA\GRL_GetCERES_year2010_2014.mat');
X = log10(DAT.lam_interp/1e3); Xmax = max(X(:));
Y = log10(nanmean(DAT.Pow.mse(:,:,TT),3).*...
    nanmean((Lv^(-2))*DAT.VAR.mse(:,TT),2)./...
    (DAT.lam_interp'.*trapz(1./DAT.lam_interp,nanmean(DAT.Pow.mse(:,:,TT),3)))); hold on;
P3 = plot(X,Y,'Linewidth',lw,'color',[0.5 0.5 1],'Linestyle','-'); hold on;

if av6==0
% Plot ERA5 reanalysis
load('MAT_DATA/UCP4_2252019_year2010_2014_month1_12.mat','DAT');
X = log10(DAT.lam_interp/1e3); Xmax = max(X(:));
Y = log10(nanmean(DAT.Pow.mse(:,:,TT),3).*...
    nanmean((Lv^(-2))*DAT.VAR.mse(:,TT),2)./...
    (DAT.lam_interp'.*trapz(1./DAT.lam_interp,nanmean(DAT.Pow.mse(:,:,TT),3)))); hold on;
P4 = plot(X,Y,'Linewidth',lw,'color','k'); hold on;
else
    % Plot ERA5 reanalysis where fields are averaged over 1 day first
load(['MAT_DATA\Test3252019_ERA_2010_2014_month1_12av6_',num2str(av6),'.mat'],'DAT');
X = log10(DAT.lam_interp/1e3);
Y = log10(nanmean(DAT.Pow.mse(:,:,TT),3).*...
    nanmean((Lv^(-2))*DAT.VAR.mse(:,TT),2)./...
    (DAT.lam_interp'.*trapz(1./DAT.lam_interp,nanmean(DAT.Pow.mse(:,:,TT),3)))); hold on;
P4 = plot(X,Y,'Linewidth',lw,'color','k'); hold on;
end

xlim([Xmin Xmax]);
set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex');
ylabel('$\varphi_{H}/\lambda\ \left[\mathrm{kg^{2}\ m^{-4}}\right]$','Fontsize',fz,'Interpreter','Latex');
xlabel('$\lambda\ \left[\mathrm{km}\right]$','Fontsize',fz,'Interpreter','Latex');
title('$\mathrm{\left(b\right)\ \lambda^{-1}\times MSE\ Power\ Spectrum\ for\ all\ datasets}$','Fontsize',fz,'Interpreter','Latex');
LEG2 = legend([P1(1) P1(2) P1(3) P2(1) P2(2) P2(3) P3 P4],...
    {'LC UNI-RAD','LC UNI-SEF','LC CTRL','NG UNI-RAD','NG UNI-SEF','NG CTRL','CERES','ERA'},...
    'Location','Eastoutside','Interpreter','Latex','Fontsize',fz);


S(1).XLim = [Xmin Xmax];

marl = 0.1; marr = 2.5e-2; marv = 5e-2; leg = 0.225;
S(1).Position = [marl 0.5+marv 1-marl-marr-leg 0.5-2*marv];
LEG1.Position = [1-marr-leg 0.5+marv leg 0.5-2*marv];
S(2).Position = [marl 2*marv 1-marl-marr-leg 0.5-2*marv];
LEG2.Position = [1-marr-leg 2*marv leg 0.5-2*marv];

G = gca; YTIK = G.YTickLabel; XTIK = G.XTickLabel;
for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
for iy = 1:numel(YTIK), YTIK{iy}=strcat('$10^{',YTIK{iy},'}$'); end
set(gca,'XTickLabel',XTIK,'YTickLabel',YTIK);

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'GRL_Fig02.pdf']); % Save plot