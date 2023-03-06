
clc 
clear
close all


%% Ch.1 Energietransitiemodel input
% originele file                merit_order.911586 t/m 75 is output; Annual load: 
% nieuwe exoprt 2 Jan 2023:     merit_order.984290 t/m 92 is output; Annual load: 184.6 TWh/year


% import 'NL Power usage (Load)' as an CSV file:
merit_order_ETM_rawimport = readtable('merit_order.984290.csv'); % ETM model - column 2 till 76 is .output - column 77 to 182 is .input powers - all in MW
% converter to matrix for datetime and matrix for doubles
time_array_orig = merit_order_ETM_rawimport{:,1};
time_array = time_array_orig + 365*2 +1;
time_array(:,2) = time_array_orig + 365*10+3;
merit_order_ETM_data = merit_order_ETM_rawimport{:,2:end};
% split producer from consumer:
m_o_producer = merit_order_ETM_rawimport(:,2:92);   % dit zijn .output
m_o_consumer = merit_order_ETM_rawimport(:,93:end); % dit zijn .input


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
Produce_curve = sum(m_o_producer{:,2:end},2); % [MW] data per hour


%% Construct consume curve
Consume_curve_orig = sum(m_o_consumer{:,2:end},2);

Cons_annual_source_TWh = sum(Consume_curve_orig)/1e6 % ETM source for 2019 is: 124.4 TWh/year which is a lot actually.
