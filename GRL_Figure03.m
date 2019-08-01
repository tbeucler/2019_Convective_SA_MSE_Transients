%% GRL_Figure03.m
% tbeucler - 6/18/2019
% Edited for Github repository of GRL submission
% tbeucler - 3/17/2019

close all; fclose('all'); clearvars;

%% 1. Parameters and initialization

% Fields and titles in each subplot
f = {'lw','sw','sef','adv'}; % Field name in structure
LEG = {'LC UNI-RAD','LC UNI-SEF','LC CTRL','NG UNI-RAD',...
    'NG UNI-SEF','NG CTRL','CERES','ERA'}; % Legend string
TIT = {'(a) Longwave Radiation','(b) Shortwave Radiation',...
    '(c) Surface Enthalpy Fluxes','(d) MSE Advection'}; % Title
XLAB = '$\lambda\ \left[\mathrm{km}\right]$'; % Units of x-axis
YLAB = '$\left[1/\mathrm{day}\right]$'; % Units of y-axis

% General appearance
facNONAGG = 5; % Divide non-aggregated simulations by constant factor
fz = 14; % Fontsize
lw = 3; % Linewidth
marg = 0.08; % Left margin 
margmid = 2.5e-2; % Middle margin
margv = 0.05; % Vertical margin
wid = 0.5-marg-margmid; % Subplot's width
widv = 0.5-2*margv; % Subplot's height

% Appearance of each line
col = [[237 177 33];[237 177 33];[237 177 33];[0 127.5 0];[0 127.5 0];...
    [0 127.5 0];[127.5 127.5 255];[0 0 0]]/255; % Color of line
ls = {':','--','-',':','--','-','-','-'}; Nsim = length(ls); % Linestyle
sim = {'LC_UNI_RAD_1day_av','LC_UNI_SEF_1day_av','LC_CTRL_1day_av',...
    'NG_UNI_RAD','NG_UNI_SEF','NG_CTRL','CERES','ERA5_1day_av'};
tmin = [40 40 40 30 30 30 0 0];
tmax = [80 80 80 1e4 1e4 1e4 1e4 1e4];

% Axes bounds
Xmin = 2; % log10(min(wavelength[km]))
Xmax = 4.2992; % log10(max(wavelength[km]))

% Physical constants
spd = 24*3600; % Number of seconds in one day

%% 2. Figure
figure
set(gcf,'Position',[50 50 1000 700]);

for isim = 1:Nsim, load(['MAT_DATA',filesep,sim{isim}]); % Load data
    N = numel(DAT.lam_interp); % Resolution of wavelength vector
    for isub = 1:4, S(isub) = subplot(2,2,isub);
        if isim==1 % Plot the zero line for reference
            line([Xmin Xmax],[0 0],'color',[0.6 0.6 0.6],'Linewidth',lw/2);
            hold on;
        end; t = DAT.t - DAT.t(1); % Time coordinate
        % Select time period to average over
        [~,i1] = min(abs(t-tmin(isim))); [~,i2] = min(abs(t-tmax(isim)));
        TT = i1:i2; SIZ = size(DAT.Agg.mse); % Time period and array's size
        % Calculate MSE power spectram and diabatic spectral rate
        if SIZ(1)==1, PHI = permute(DAT.Agg.mse,[2 1 3])/2;
            NUM = permute(DAT.Agg.(f{isub})(:,:,TT),[2 1 3]);
        else, PHI = DAT.Agg.mse/2; NUM = DAT.Agg.(f{isub})(:,:,TT);
        end
        % Find normalization constant int{phi d(1/lambda)})
        Iphi = trapz(1./DAT.lam_interp,PHI,1);
        % Tranpose the wavelength vector for the x-axis if necessary
        SIZ = size(DAT.lam_interp);
        if SIZ(1)>1, DAT.lam_interp = DAT.lam_interp'; end
        % Use logarithm 10 of wavelength [km] for x-axis
        X = log10(DAT.lam_interp/1e3);
        % Use normalized tendency of variance [1/day] for y-axis
        Y = nanmean(spd*NUM./(repmat(DAT.lam_interp',1,1,numel(TT)).*...
            repmat(Iphi(1,1,TT),N,1,1)),3); hold on;
        % Divide the spectral rate by a factor facNONAGG if non-aggregated
        if isim==1||isim==4, Y = Y/facNONAGG; end
        % Plot using custom color anda linestyle
        PLOT(isim) = plot(X,Y,'Linewidth',lw,'color',col(isim,:),...
            'Linestyle',ls{isim}); hold on;
        % Subplot's appearance
        xlim([Xmin Xmax]); % x bounds
        ylabel(YLAB,'Fontsize',fz,'Interpreter','Latex');
        title(TIT{isub},'Fontsize',fz,'Interpreter','Latex');
        set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex','Xdir','reverse');
        % Adjust location of y-axis
        if isub==2||isub==4, set(gca,'YAxisLocation','right'); end
        % Create legend
        if isim==Nsim, L = legend(PLOT,LEG,'Location','EastOutside',...
                'Orientation','vertical','Interpreter','Latex');
        end
        % Logarithmic scale for x-axis
        if isub>=3 && isim==4, ax = gca;
            G = gca; XTIK = G.XTickLabel;
            for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
            set(ax,'XTickLabel',XTIK); set(ax.XLabel,...
                'Interpreter','Latex','String',XLAB,'Fontsize',fz);
        elseif isub<3, set(gca,'XTickLabel','');
        end
    end
end

% Adjust subplot's position
S(1).Position = [marg 0.5+margv wid widv];
S(2).Position = [0.5+margmid 0.5+margv wid widv];
S(3).Position = [marg 1.5*margv wid widv];
S(4).Position = [0.5+margmid 1.5*margv wid widv];
% Adjust legend's position
set(L,'Position',[S(1).Position(1)+0.6*S(1).Position(3) ...
    S(1).Position(2)+0.4*S(1).Position(4) 0.4*S(1).Position(3) ...
    0.4*S(1).Position(4)]);

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'GRL_Fig03.pdf']); % Save plot