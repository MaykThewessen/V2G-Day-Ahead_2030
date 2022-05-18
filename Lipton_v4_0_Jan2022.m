% Origin:       Tuesday 8 Feb 2022
% Author:       Mayk Thewessen
% Department:   Strategy - Research
% Intent:       Electricity market NL analysis for Vehicle-to-Grid in 2030


close all
clear
clc

format short eng
set(groot,'defaultLineLineWidth',2)


% %% import 'NL Power usage (Load)' as an CSV file:
% M = readtable('merit_order.891876.csv');
% 
% %% import xlsx file 
% Load_xlsx_readtable = readtable('Total Load - Day Ahead _ Actual_202101010000-202201010000.xlsx');
% 
% %% Import Electrical Load of NL in 2021 whole year - actual - XML
% Load_readtable_1year = readtable('ACTUAL_TOTAL_LOAD_202101010000-202201010000.xml');
% plot(Load_readtable_1year.quantity)
% 
% plot(Load_readtable_1year.quantity(:))
% 
% %%
% A_column_vector = (Load_readtable_1year.quantity(:));
% A_row_vector = (Load_readtable_1year.quantity(:)');
% 
% %% option 2
% Load_readtable_1day = readtable('ACTUAL_TOTAL_LOAD_202102040000-202102050000.xml');


%% Ch.1 Energietransitiemodel input

% import 'NL Power usage (Load)' as an CSV file:
merit_order_ETM_rawimport = readtable('merit_order.904566.csv'); % ETM model - column 2 till 76 is .output - column 77 to 182 is .input powers - all in MW
% converter to matrix for datetime and matrix for doubles
time_array = merit_order_ETM_rawimport{:,1};
merit_order_ETM_data = merit_order_ETM_rawimport{:,2:end};
% split producer from consumer:
m_o_producer = merit_order_ETM_rawimport(:,2:75);
m_o_consumer = merit_order_ETM_rawimport(:,76:end);


% Find largest contributors to producers
merit_prod_summed_per_type = sum(m_o_producer{:,:});
[merit_prod_max, merit_prod_pos] = sort(merit_prod_summed_per_type,'descend');                 % prod max in MWh per year, position is 1 column number less than in rawimport
find_top = 10; % generators that contribute yearly most volume in MWh/year
merit_prod_pos_top = merit_prod_pos(:,1:find_top)+1;
merit_prod_max_top = merit_prod_max(:,1:find_top);
% A = convertCharsToStrings(merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top+1})


% Find largest contributors to producers
merit_cons_summed_per_type = sum(m_o_consumer{:,:});
[merit_cons_max, merit_cons_pos] = sort(merit_cons_summed_per_type,'descend');                 % prod max in MWh per year, position is 1 column number less than in rawimport
find_top = 10; % generators that contribute yearly most volume in MWh/year
merit_cons_pos_top = merit_cons_pos(:,1:find_top)+1+74;
merit_cons_max_top = merit_cons_max(:,1:find_top);





%% Construct demand curve 2030
Produce_curve = sum(m_o_producer{:,2:end},2);



%% Construct consume curve
Consume_curve = sum(m_o_consumer{:,2:end},2);




%% Step 2: Calculate price
% 2 = PV buildings rooftop solar (industry)
% 44 = PV large scale solar
% 57 = coastal wind energy
% 58 = inland wind energy
% 59 = offshore wind energy
% 61 = PV households
% A) subtract no marginal cost from power usage curve: construct residual load curve
PV_sum_producers = merit_order_ETM_rawimport{:,2}+merit_order_ETM_rawimport{:,44}+merit_order_ETM_rawimport{:,61} ; 
Wind_sum_producers = merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59} ;
renewable_producers = merit_order_ETM_rawimport{:,2}+merit_order_ETM_rawimport{:,44}+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59}+merit_order_ETM_rawimport{:,61};
residual_load_curve = Consume_curve - renewable_producers;
price_electricity = 21.486.*exp(residual_load_curve.*1e-4) - 21.486; % v2: y = 21,486e0,0001x, v1: y = 27,775e4E-05x



%% initiliaze plot
h0 = figure('Name','Electricity market NL 2030','pos',[0 0 2000 1200]); % width and height start and end points
subplot_size_x = 2; % height - number of rows
subplot_size_y = 2; % width - number of columns
subplot_count = 0;

%% Plot A
subplot(subplot_size_x,subplot_size_y,    0 * subplot_size_y + 1) % x*sub_y bepaalt row, sub_y +x bepaalt welke column
plot(time_array,Produce_curve)
hold on
plot(time_array, Consume_curve)
plot(time_array,residual_load_curve) % Why is the residual load sometimes negative? thus more renwable producers than consumers? or renewables that do negative power?
xlabel('Time')
ylabel('Electrical Power [MW]')
grid
title('Whole year')

yyaxis right
plot(time_array,price_electricity)
legend('Production total','Consumption total','Residual load','Electricity price')
ylabel('€/MWh')

%% Plot B
subplot(subplot_size_x,subplot_size_y, 1 * subplot_size_y + 1)
% sort electricity price
duration_curve_price = sort(price_electricity,'descend');
area(1:8760,duration_curve_price)

%histogram(price_electricity)
ylabel('Electricity price €/MWh')
xlabel('Occurance [hours per year]')
grid
% set(gca,'YScale','log')


%% Duration curve power
duration_curve_power = sort(Consume_curve,'descend')/1000;
%duration_curve_power((length(duration_curve_power)+1),1) = 0;
area(1:length(duration_curve_power),duration_curve_power)
ylim([0 max(duration_curve_power)])
xlabel('Hours per year')
ylabel('Electrical power [GW]')
grid
legend('Duration curve consumption')
title('Whole year')



%% Plot C rechtsboven - Production and consumption of one week - wind energy power into graph
subplot(subplot_size_x,subplot_size_y, 1 * subplot_size_y + 0)
plot(time_array,Produce_curve)
hold on
plot(time_array, Consume_curve)
xlabel('Time')
ylabel('Electrical Power [MW]')
title('Summer week')

% choose time
start_point = 2500; % 7 June 2030
days = 7;
xlim([time_array(start_point), time_array(start_point+days*24)])

% Solar (everything):
area(time_array, residual_load_curve+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59}+merit_order_ETM_rawimport{:,44}+merit_order_ETM_rawimport{:,2}+merit_order_ETM_rawimport{:,61}) % wind + stacked: large scale solar

% Fossil as second foreground
area(time_array, Consume_curve - PV_sum_producers)
% area(time_array, residual_load_curve+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59})

% Wind: as last foreground color
area(time_array, Wind_sum_producers)
legend('Production total','Consumption total','PV solar (household+buildings+central)','Residual load (mainly fossil backup)','Wind energy (inland, coastal, and offshore)')
grid

% ik wil graag overshot ook laten zien met stippelijn erboven over, of negatief?
% negatief: opladen van batterij


%%
% Figure
% area(time_table, Consume_curve)
% area(time_table, merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59}+merit_order_ETM_rawimport{:,44}+merit_order_ETM_rawimport{:,2}+merit_order_ETM_rawimport{:,61}) % wind + stacked: large scale solar
% area(time_table, merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59}) % all wind energy
% area(residual_load_curve)
% %plot(time_table,residual_load_curve)
% legend('Production total','Consumption total','PV solar (household+buildings+central)','Wind energy (inland, coastal, and offshore)','Residual load curve')
% grid

% Het is eigenlijk mooier om PV en Wind tegen bovenkant van consumption curve aan te plakken


%% Plot D - Residual load curve and Price curve
subplot(subplot_size_x,subplot_size_y, 1 * subplot_size_y + 2)
plot(time_array,Consume_curve)
hold on
plot(time_array,residual_load_curve)
xlim([time_array(start_point), time_array(start_point+days*24)])
grid
xlabel('Time')
ylabel('Electrical Power [MW]')

yyaxis right
plot(time_array,price_electricity)
hold on
xlim([time_array(start_point), time_array(start_point+days*24)])
legend('Consume curve','Residual load','Electricity price')
ylabel('€/MWh')
title('Electricity prices')


%% Find names of largest producers
for j = 1:find_top
    Producer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top(j)}));   % first column: name
    Consumer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_cons_pos_top(j)}));   % first column: name
end
Produced_MWh_year = merit_prod_max_top';
Consumed_MWh_year = merit_cons_max_top';
Largest_prod_cons_table = table(Producer_type_name,Produced_MWh_year,merit_prod_pos_top',Consumer_type_name,Consumed_MWh_year,merit_cons_pos_top');


%% import power plant capacity values
power_plant_capacities = readtable('production_parameters.895094.csv');









