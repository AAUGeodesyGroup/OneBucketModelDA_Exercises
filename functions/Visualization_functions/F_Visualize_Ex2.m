function F_Visualize_Ex2(P_ens,S_ens,Q_ens,sim)

%% Visualization

clr_ens  = [0.5 0.5 0.5]; % color ensemble
lw = 2;
msize = 15;
fsize = 12;
opt_save = 1;
Ne = size(P_ens,2);

% precipitation
figure
set(gca,'FontSize',fsize)
hold on
for ii=1:size(P_ens,2)
    h1 = plot(sim(1):sim(2),P_ens(sim(1):sim(2),ii),'.-','Color',clr_ens,'LineWidth',lw,'MarkerSize',msize);
end
h2 = plot(sim(1):sim(2),mean(P_ens(sim(1):sim(2),:),2),'.-','Color','k','LineWidth',lw,'MarkerSize',msize);
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('P-E'); 
title('Perturbed input data: precipitation - evapotranspiration ensemble');
legend([h2 h1],'Ensemble mean','Ensemble members')
box on
if(opt_save==1) print(gcf,'figures/Ex2_P_OL','-djpeg'); end
% storage
figure
set(gca,'FontSize',fsize)
hold on
for ii = 1:Ne
    h1 = plot(sim(1):sim(2),S_ens(:,ii),'.-','Color',clr_ens,'LineWidth',lw,'MarkerSize',msize);
end
h2 = plot(sim(1):sim(2),mean(S_ens,2),'.-','Color','k','LineWidth',lw,'MarkerSize',msize);
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('S')
title('Water storage estimates - model simulation ensemble');
legend([h2 h1],'Ensemble mean','Ensemble members')
box on
ylim([0 max(max(S_ens))+2])
if(opt_save==1) print(gcf,'figures/Ex2_S_OL','-djpeg'); end
% discharge
figure
set(gca,'FontSize',fsize)
hold on
for ii = 1:Ne
    h1 = plot(sim(1):sim(2),Q_ens(:,ii),'.-','Color',clr_ens,'LineWidth',lw,'MarkerSize',msize);
end
h2 = plot(sim(1):sim(2),mean(Q_ens,2),'.-','Color','k','LineWidth',lw,'MarkerSize',msize);
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('R')
title('River discharge estimates - model simulation ensemble');
box on
ylim([0 max(max(S_ens))+2])
legend([h2 h1],'Ensemble mean','Ensemble members')
if(opt_save==1) print(gcf,'figures/Ex2_R_OL','-djpeg'); end

end

