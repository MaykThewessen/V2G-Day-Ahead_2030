close all
clc
clear

load('ETM_3_March')

Vehicles = ETM_3_March(:,1);
Revenue = ETM_3_March(:,2);
Costs = ETM_3_March(:,3);

V2G_overview = table(Vehicles,Revenue,Costs)

Vehicles_million = Vehicles./1e6;

%% Plot
h0 = figure('pos',[1200 900 400 250]);

bar(Vehicles_million,[Revenue, Costs], 'BarWidth', 2)

grid


xticks(Vehicles_million)  % xticks([round(Vehicles_million,1)])
xtickformat('%.1f M')
ytickformat('€%.0f')

%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%xticks([linspace(0,3e6,7)])
%xticks(Vehicles);
%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%xtickformat('%g')


legend('Electricity revenue','Electricity costs')

xlabel('Number of vehicles connected to V2G pole [million]')
ylabel('Annual average V2G Electricity price [€/MWh]')

ylim([0 max(Revenue)*1.15])

save_fig(h0,'ETM_revenue_costs_plot_v2');
