%% GRL_Figure03.m
% tbeucler - 3/17/2019

close all; fclose('all'); clearvars;

COL = [[7 118 190];[218 89 33];[237 177 33];[127 48 143];[127 177 61]]/255; % MSE power spectra at different times
cmap = [117 112 179; 217 150 0; 231 41 138; 27 158 119]/255; % LW,SW,SEF,ADV colors
green = [0 0.5 0];
fz = 14; lw = 3; % Fontsize and Linewidth
Xmin = 5-3; Xmax = 7.2992-3;

av6 = 24;
facNONAGG = 5; % Divide non aggregated simulations by constant factor
N = 100;
spd = 24*3600;

f = {'lw','sw','sef','adv'};
TIT = {'(a) Longwave Radiation','(b) Shortwave Radiation','(c) Surface Enthalpy Fluxes','(d) MSE Advection'};

figure
set(gcf,'Position',[50 50 1000 700]);
Nsim = 8;
for isim = 1:Nsim
    if isim==1, SST = 1; rad = 'cam';
        if av6>0, filename = ['UCP2_2212109_SST',num2str(SST),'_rad',rad,'_av6',num2str(av6)];
            Ladv = load(['MAT_DATA/UCP2_2212109_SST',num2str(SST),'_rad',rad],'DAT');
            DAT.Agg.adv = Ladv.DAT.Agg.adv;
        else, filename = ['UCP2_2212109_SST',num2str(SST),'_rad',rad];
        end; load(['MAT_DATA/',filename],'DAT');
        col = COL(3,:); ls = ':'; t1 = 50; t2 = 80;
        DAT.Agg.lw = 0*DAT.Agg.lw.^0; DAT.Agg.sw = 0*DAT.Agg.sw.^0;
        DAT.Agg.sef = DAT.Agg.sef/facNONAGG;
        DAT.Agg.adv = DAT.Agg.adv/facNONAGG;
        tend_lw = 0;
    elseif isim==2, SST = 2; rad = 'cam';
        if av6>0, filename = ['UCP2_2212109_SST',num2str(SST),'_rad',rad,'_av6',num2str(av6)];
            Ladv = load(['MAT_DATA/UCP2_2212109_SST',num2str(SST),'_rad',rad],'DAT');
            DAT.Agg.adv = Ladv.DAT.Agg.adv;
        else, filename = ['UCP2_2212109_SST',num2str(SST),'_rad',rad];
        end; load(['MAT_DATA/',filename],'DAT');
        col = COL(3,:); ls = '--';
        DAT.Agg.sef = 0*DAT.Agg.sef.^0;
    elseif isim==3, SST = 300; rad = 'cam';
        if av6>0, filename = ['UCP2_2212109_SST',num2str(SST),'_rad',rad,'_av6',num2str(av6)];
            Ladv = load(['MAT_DATA/UCP2_2212109_SST',num2str(SST),'_rad',rad],'DAT');
            DAT.Agg.adv = Ladv.DAT.Agg.adv;
        else, filename = ['UCP2_2212109_SST',num2str(SST),'_rad',rad];
        end; load(['MAT_DATA/',filename],'DAT');
        col = COL(3,:); ls = '-';
        tend_lw = trapz(1./DAT.lam_interp,nanmean(DAT.Agg.dmsedt,3))/...
            trapz(1./DAT.lam_interp,nanmean(DAT.Agg.lw,3));
    elseif isim==4, load('MAT_DATA/GRL_GetNG_2','DAT');
        col = green; ls = ':'; t1 = 30; t2 = 99999;
        DAT.Agg.lw = 0*DAT.Agg.lw.^0; DAT.Agg.sw = 0*DAT.Agg.sw.^0;
        DAT.Agg.sef = DAT.Agg.sef/facNONAGG;
        DAT.Agg.adv = DAT.Agg.adv/facNONAGG;
    elseif isim==5, load('MAT_DATA/GRL_GetNG_3','DAT');
        col = green; ls = '--';
        DAT.Agg.sef = 0*DAT.Agg.sef.^0;
    elseif isim==6, load('MAT_DATA/GRL_GetNG_1','DAT');
        col = green; ls = '-';
        tend_lw = trapz(1./DAT.lam_interp,nanmean(DAT.Agg.dmsedt,3))/...
            trapz(1./DAT.lam_interp,nanmean(DAT.Agg.lw,3));
    elseif isim==8
        if av6>0, load(['MAT_DATA\Test3252019_ERA_2010_2014_month1_12av6_',num2str(av6),'.mat'],'DAT');
            Ladv = load('MAT_DATA/UCP4_2252019_year2010_2014_month1_12.mat','DAT');
            DAT.Agg.adv = Ladv.DAT.Agg.adv;
        else, load('MAT_DATA/UCP4_2252019_year2010_2014_month1_12.mat','DAT');
        end
        col = 'k'; ls = '-'; t1 = 0; t2 = 9999;
        tend_lw = trapz(1./DAT.lam_interp,nanmean(DAT.Agg.dmsedt,3))/...
            trapz(1./DAT.lam_interp,nanmean(DAT.Agg.lw,3));
    elseif isim==7, load('MAT_DATA\GRL_GetCERES_year2010_2014.mat');
        col = [0.5 0.5 1]; ls = '-'; t1 = 0; t2 = 9999;
        DAT.Agg.sef = NaN*DAT.Agg.lw.^0;
        DAT.Agg.adv = NaN*DAT.Agg.lw.^0;
    elseif isim==9, load('MAT_DATA\Test3252019_ERA_2010_2014_month1_12av6_24.mat','DAT');
        col = [0 0 0.75]; ls = '-'; t1 = 0; t2 = 9999;
        DAT.Agg.adv = NaN*DAT.Agg.lw.^0;
    end
    disp(['Time-averaged MSE tendency divided by lw is ',num2str(tend_lw)]);
    
    for isub = 1:4, S(isub) = subplot(2,2,isub);
        if isim==1, line([Xmin Xmax],[0 0],'color',[0.6 0.6 0.6],'Linewidth',lw/2); hold on; end
        
        t = DAT.t - DAT.t(1);
        [~,i1] = min(abs(t-t1)); [~,i2] = min(abs(t-t2)); TT = i1:i2;
        if isim>3&&isim<7, PHI = permute(DAT.Agg.mse,[2 1 3])/2;
            NUM = permute(DAT.Agg.(f{isub})(:,:,TT),[2 1 3]);
        else, PHI = DAT.Agg.mse/2; NUM = DAT.Agg.(f{isub})(:,:,TT);
        end
        Iphi = trapz(1./DAT.lam_interp,PHI,1); % Find normalization constant int{phi d(1/lambda)})
        
        X = log10(DAT.lam_interp/1e3);
        Y = nanmean(spd*NUM./(repmat(DAT.lam_interp',1,1,numel(TT)).*...
            repmat(Iphi(1,1,TT),N,1,1)),3); hold on;
        
        PLOT(isim) = plot(X,Y,'Linewidth',lw,'color',col,'Linestyle',ls); hold on;
        xlim([Xmin Xmax]);
        ylabel('$\left[1/\mathrm{day}\right]$','Fontsize',fz,'Interpreter','Latex');
        title(TIT{isub},'Fontsize',fz,'Interpreter','Latex');
        set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex','Xdir','reverse');
        if isub==2||isub==4, set(gca,'YAxisLocation','right'); end
        if isim==Nsim
            L = legend(PLOT,{'LC UNI-RAD','LC UNI-SEF','LC CTRL','NG UNI-RAD','NG UNI-SEF','NG CTRL','CERES','ERA'},...
                'Location','EastOutside','Orientation','vertical','Interpreter','Latex');
        end
        if isub>=3 && isim==4, ax = gca;
            G = gca; XTIK = G.XTickLabel;
            for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
            set(ax,'XTickLabel',XTIK);
            ax.XLabel.Interpreter = 'Latex';
            ax.XLabel.String = '$\lambda\ \left[\mathrm{km}\right]$';
            ax.XLabel.FontSize = fz;
        elseif isub<3, set(gca,'XTickLabel','');
        end
    end
end

marg = 0.08; margmid = 2.5e-2; margv = 0.05; wid = 0.5-marg-margmid; widv = 0.5-2*margv;
S(1).Position = [marg 0.5+margv wid widv];
S(2).Position = [0.5+margmid 0.5+margv wid widv];
S(3).Position = [marg 1.5*margv wid widv];
S(4).Position = [0.5+margmid 1.5*margv wid widv];
set(L,'Position',[S(1).Position(1)+0.6*S(1).Position(3) S(1).Position(2)+0.4*S(1).Position(4) 0.4*S(1).Position(3) 0.4*S(1).Position(4)]);

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'GRL_Fig03.pdf']); % Save plot