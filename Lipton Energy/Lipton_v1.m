% Origin:       Tuesday 8 Feb 2022
% Author:       Mayk Thewessen
% Department:   Strategy - Research
% Intent:       Electricity market NL analysis for Vehicle-to-Grid in 2030
% 
% Show steady state TMS+HVAC consumption relative to WLTP usage numbers
% in Wh/km. Goal: to show the influence of ambient temperature on usage and
% range, and to show how much impact can be made.
%
% part 1: show steady state TMS usage vs 85 Wh/km
% part 2: add dynamic heat up, cool down effect
% part 3: add 3 climates; Trondheim, Amsterdam, Cordoba; show histogram
% part 4: show speed dependencies, 46.5 km/h vs 100 km/h vs 130 km/h
% part 5: improve with accurate vehicle usage vs amb temp dependency of air
% density and thus air resistance and by rolling resistance temperature influence.


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
merit_order_ETM_rawimport = readtable('merit_order.895094.csv'); % ETM model - column 2 till 76 is .output - column 77 to 182 is .input powers - all in MW
% converter to matrix for datetime and matrix for doubles
merit_order_ETM_time_table = merit_order_ETM_rawimport{:,1};
merit_order_ETM_data = merit_order_ETM_rawimport{:,2:end};
% split producer from consumer:
m_o_producer = merit_order_ETM_rawimport(:,2:end);
m_o_consumer = merit_order_ETM_rawimport(:,76:end);


% Find largest contributors to producers
merit_prod_summed = sum(m_o_producer{:,:},2);
[merit_prod_max, merit_prod_pos] = sort(merit_prod_summed,'descend');                 % prod max in MWh per year, position is 1 column number less than in rawimport
find_top = 10; % generators that contribute yearly most volume in MWh/year
merit_prod_pos_top = merit_prod_pos(:,1:find_top);
merit_prod_max_top = merit_prod_max(:,1:find_top);
% A = convertCharsToStrings(merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top+1})


% Find largest contributors to producers
merit_cons_summed = sum(m_o_consumer{:,:},2);
[merit_cons_max, merit_cons_pos] = sort(merit_cons_summed,'descend');                 % prod max in MWh per year, position is 1 column number less than in rawimport
find_top = 10; % generators that contribute yearly most volume in MWh/year
merit_cons_pos_top = merit_cons_pos(:,1:find_top);
merit_cons_max_top = merit_cons_max(:,1:find_top);





%% Construct demand curve 2030
Produce = sum(m_o_producer(:,2:end),2);



%% Construct consume curve
Consume = sum(m_o_consumer(:,2:end),2);





%% Find names of largest producers
for j = 1:find_top
    Producer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top(j)+1}));   % first column: name
    Consumer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_cons_pos_top(j)+1}));   % first column: name
end
Produced_MWh_year = merit_prod_max_top';
Consumed_MWh_year = merit_cons_max_top';
Largest_prod_cons_table = table(Producer_type_name,Produced_MWh_year,Consumer_type_name,Consumed_MWh_year);


%% construct Load = added up all generation and flex
% Load_ETM = sum(merit_order_ETM_data(:,2:end),2);


%% import power plant capacity values
power_plant_capacities = readtable('production_parameters.895094.csv');









