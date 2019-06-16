%% GRL_Figure04.m
% tbeucler - 3/17/2019

close all; fclose('all'); clearvars;

fz = 16; lw = 3;

Lv = 2.5e6;

LCSSTa = [1 2 300];
INV = [1 3 2 4]; % Inverts panels (b) and (c)
leva = [850 750 500 250]; Nlev = numel(leva);
LEG = {'850hPa','750hPa','500hPa','250hPa'};
TIT = {'(a) Spectra at Fixed Vertical Level','(c) Ratio of Coherences','(b) Self-Aggregation Number','(d) Ratio of Variances'};

Xmin = 3.7782; Xmax = 7.0895;

addpath(genpath('cbrewer'));
COL = cbrewer('qual','Set1',2*Nlev);
COL(3,:) = COL(1,:);
COL(1,:) = [237 177 33]/255; % Same yellow as LC

figure
set(gcf,'Position',[50 50 1000 750]);

for isim = [1 3], LCSST = LCSSTa(isim);
    load(['C:\Users\Tom\Desktop\Radiative_convective_instability\Cloud_rad_water_vap\MAT_DATA\GRLFig4LC_SST',...
        num2str(LCSST),'_tini360000_Nstep120.mat']);
    if isim==1, ls = ':';
    elseif isim==2, ls = '--';
    elseif isim==3, ls = '-';
    end
    for isub = 1:4, S(INV(isub)) = subplot(2,2,isub);
        if isim==3&&isub~=4
            line([Xmin Xmax],[0 0],'color',[0.6 0.6 0.6],'Linewidth',lw/2); hold on;
        end
        
        for ilev = 1:Nlev, LEV = ['lev',num2str(leva(ilev))];
            
            if isub==1, vari = log10(nanmean(DAT.(LEV).Pow.qv,3).*...
                    nanmean(DAT.(LEV).VAR.qv/1e6,3)./...
                    (DAT.(LEV).lamoct.*...
                    trapz(1./DAT.(LEV).lamoct,nanmean(DAT.(LEV).Pow.qv,3))));
                var1 = vari;
                var2 = log10(nanmean(DAT.(LEV).Pow.dse,3).*...
                    nanmean(DAT.(LEV).VAR.dse/Lv^2,3)./...
                    (DAT.(LEV).lamoct.*...
                    trapz(1./DAT.(LEV).lamoct,nanmean(DAT.(LEV).Pow.dse,3))));
                YLIM = [-9.9 -5.9];
            elseif isub==3, var1 = nanmean(DAT.(LEV).FIG.wradq,3)./...
                    nanmean(DAT.(LEV).FIG.wnonq,3); var2 = var1; vari = var1;
                YLIM = [-0.29 0.29];
            elseif isub==2, var1 = nanmean(DAT.(LEV).FIG.cohwradq,3);
                var2 = nanmean(DAT.(LEV).FIG.cohwnonq,3);
                vari = var1./var2; YLIM = [-3.9 7.9];
            elseif isub==4, var1 = log10(nanmean(DAT.(LEV).FIG.wrad,3));
                var2 = log10(nanmean(DAT.(LEV).FIG.wnon,3));
                vari = var1./var2; YLIM = [-0.7 0.7];
            end
            
            X = log10(DAT.(LEV).lamoct);
            Y1 = var1; Y2 = var2; Y = vari;
            %plot(X,Y1,'color',COL(ilev,:),'Linewidth',lw,'LineStyle',ls); hold on;
            %plot(X,Y2,'color',COL(ilev,:),'Linewidth',lw,'LineStyle',ls); hold on;
            P(isim,ilev) = plot(X,Y,'color',COL(ilev,:),'Linewidth',lw,'LineStyle',ls); hold on;
            if isub==1&&isim==3, Q(ilev) = plot(X,Y2,'color',COL(ilev,:),'Linewidth',lw,'LineStyle','--'); hold on; end
            title(TIT(isub),'Fontsize',fz,'Interpreter','Latex');
            xlim([Xmin Xmax]);
            ylim(YLIM);
            
            
            
        end
        set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex','Xdir','reverse'); 
        if isub==1||isub==3, grid on; end
        if isub==3&&isim==3, L2 = legend(P(3,:),LEG,'Location','Southoutside','Orientation','vertical',...
                'Fontsize',fz,'Interpreter','Latex'); 
        elseif isub==1&&isim==3, L = legend([P(3,1) P(1,1) Q(1)],...
                {'$\varphi_{q}\ \mathrm{CTRL}$','$\varphi_{q}\ \mathrm{UNIRAD}$','$\varphi_{s}\ \mathrm{CTRL}$'},...
                'Orientation','vertical','Fontsize',fz,'Interpreter','Latex');
        end
    end
end

marg = 0.08; margmid = 2.5e-2; margv = 0.05; wid = 0.5-marg-margmid; widv = 0.5-2*margv;
S(1).Position = [marg 0.5+margv wid widv];
S(2).Position = [0.5+margmid 0.5+margv wid widv];
S(3).Position = [marg 1.5*margv wid widv];
S(4).Position = [0.5+margmid 1.5*margv wid widv];
set(L2,'Position',[S(3).Position(1)+0.6*S(3).Position(3) S(3).Position(2)+0.6*S(3).Position(4) 0.4*S(3).Position(3) 0.4*S(3).Position(4)]);
set(L,'Position',[S(1).Position(1)+0.55*S(1).Position(3) S(1).Position(2)+0.65*S(1).Position(4) 0.45*S(1).Position(3) 0.35*S(1).Position(4)]);

for isub = 1:4, S(isub);
    
    if isub==2||isub==4, set(S(isub),'YAxisLocation','right'); end
    if isub==1
        G = S(isub); YTIK = G.YTickLabel;
        for ix = 1:numel(YTIK), YTIK{ix}=strcat('$10^{',YTIK{ix},'}$'); end
        set(S(isub),'YTickLabel',YTIK);
        G.YLabel.Interpreter = 'Latex';
    elseif isub==4
       G = S(isub); YTIK = G.YTickLabel;
        for ix = 1:numel(YTIK), YTIK{ix}=num2str(10^str2double(YTIK{ix}),'%02.2f'); end
        set(S(isub),'YTickLabel',YTIK);
        G.YLabel.Interpreter = 'Latex'; 
    end
    
    if (isub==4||isub==3) && isim==3 && ilev==Nlev, ax = S(isub);
        G = S(isub); XTIK = G.XTickLabel;
        for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
        set(ax,'XTickLabel',XTIK);
        ax.XLabel.Interpreter = 'Latex';
        ax.XLabel.String = '$\lambda\ \left[\mathrm{m}\right]$';
        ax.XLabel.FontSize = fz;
    elseif isub==1 || isub==2, set(S(isub),'XTickLabel','');
    end
    
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'GRL_Fig04.pdf']); % Save plot