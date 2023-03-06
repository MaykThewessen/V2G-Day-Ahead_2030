close all
clear
clc

set(groot,'defaultLineLineWidth',2)

%% copy ETM merit order table -> google sheets -> sort -> export as .xlsx -> import
merit_xlsx_import_generated_script % etm_merit_loaded


%% Construct data
names = etm_merit_loaded.Technology;
prices = etm_merit_loaded.MarginalCosts; % Y-axis
instal_capacity = etm_merit_loaded.InstalledCapacity./1000; % X-axis %[GW] power

instal_capacity_cumsum = cumsum(instal_capacity);


%% Plot using single bar graphs and hold function
% found idea on: https://comp.soft-sys.matlab.narkive.com/fWK4AmsP/varying-width-of-bar-graphs
instal_capacity_cumsum_from_zero = [0; instal_capacity_cumsum];

close all
h0 = figure;
hold on



newcolors = [0.9290 0.6940 0.1250 % PV
             0.9290 0.6940 0.1250
             0.9290 0.6940 0.1250
             0.9290 0.6940 0.1250

             0 0.4470 0.7410 % Wind
             0 0.4470 0.7410
             0 0.4470 0.7410

             0 0 1 % hydro

             0 0 0 % nuclear

             0.4660 0.6740 0.1880 % biomass

             0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54

          
             0.8500 0.3250 0.0980
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330
             0.6350 0.0780 0.1840
             
             
             rand rand rand
             rand rand rand
             rand rand rand
             rand rand rand
             rand rand rand
             rand rand rand
             ];




%newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
         
colororder(newcolors)

for ii = 1:length(instal_capacity_cumsum)
    center_x = (instal_capacity_cumsum_from_zero(ii) + instal_capacity_cumsum_from_zero(ii+1)) /2;
    height = prices(ii)+2;
    width = instal_capacity(ii);
    bar(center_x,height,width)
end

grid
ylabel('Marginal generation costs of electricity [€/MWh]')
xlabel('Installed generation capacity (cumulative) [GW]')
title('(at a gas TTF price of €50/MWh thermal)')

% TODO: how to make legend from a array of strings
%legend('a','b','Location','NorthWest')


legend(names,'Location','NorthWest')

save_fig(h0,'ETM_merit_order_import_make_graph_v2');






%% other idea to plot using: rectangles
% %% Plot using rectangles
% 
% % sample data: set the start of each bar, the bottom (here 0), the width and the height
% 
% x = [0; instal_capacity_cumsum(1:end-1,:)]; % start of bar
% dx = instal_capacity; % width of bar
% y = zeros(length(x),1); % start rectangle from bottom
% dy = prices; % rectanlge height stops at price
% 
% figure, hold on
% for ii=1:length(x)
%     color = [rand rand rand]; % https://www.mathworks.com/help/matlab/ref/rectangle.html
%     rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)],'FaceColor',color)
% end
% grid
% legend('a','b')
% 
% ylabel('Electricity price')
% xlabel('Generation power')
% 
% 
% %% Plot example using rectangles
% % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
% % sample data: set the start of each bar, the bottom (here 0), the width and the height
% 
% x = [0.5 0.6 0.9 1 1.2]; % start of bar
% y = zeros(length(x),1);
% dx = diff([x 1.8]); % width of bar
% dy = [1 3 2 .5 .1];
% 
% figure, hold on
% for ii=1:length(x)
%     rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)])
% end
% axis([0.5 2 0 4.1])
% 
% ylabel('Prob density')
% xlabel('Time')


%% ABC







