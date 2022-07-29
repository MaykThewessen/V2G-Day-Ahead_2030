
import_Elaad_open 


%%
h = figure('Name','Charge pole Max Power histogram','pos',[1100 800 450 300]);



h1 = histogram(MaxPower,'Normalization','probability','BinWidth',1);


xlabel('Charge power [kW]')
ylabel('Occurance')

legend('Elaad 2019 dataset of 10.000 public charges')

grid

save_fig(h,'Max_Power_hist_Elaad_v2');     % uses minimized edge borders
