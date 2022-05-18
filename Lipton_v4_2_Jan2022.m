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


%% initiliaze plot
h0 = figure('Name','Electricity market NL 2030','pos',[0 0 2000 1200]); % width and height start and end points
subplot_size_x = 4; % height - number of rows
subplot_size_y = 2; % width - number of columns
subplot_count = 0;


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
Produce_curve = sum(m_o_producer{:,2:end},2); % [MW] data per hour



%% Construct consume curve
Consume_curve = sum(m_o_consumer{:,2:end},2);


%% Import production capacities
production_parameters = readtable('production_parameters.904566.csv'); %

% Calculate installed power per type of generator:
%production_parameters.installed_power = production_parameters{:,2} .* production_parameters{:,3};
production_parameters.installed_power = production_parameters{:,"number_of_units"} .* production_parameters{:,"electricity_output_capacity_MW_"};


% Sum up the renewable contributors
P_solar = sum( production_parameters{[14,95,120],"installed_power"} );
% pv households = 120
% pv buildings = 14
% pv solar parks = 95

P_wind = sum( production_parameters{110:112,"installed_power"} );
% wind onshore 111
% wind coastal 110
% wind offshore 112

% Scale solar and wind to:
%P_zon_prognose_2030 = 33000; % [MW] laag scenario - als er veel grid congestie is - eprijs dempt flink - curtailment issues in overheidsregeling
P_zon_prognose_2030 = 46200; % [MW] hoog scenario - pv cost down
zon_scale = P_zon_prognose_2030 / P_solar

%P_wind_prognose_2030 = 8800 + 16700; % [MW] laag scenario 8.8GW onshore + 16.7GW offshore
P_wind_prognose_2030 = 8800 + 21300; % [MW] hoog scenario 8.8GW onshore + 21.3GW offshore reeds aangekodigd door overheid, plannen die dit bewerkstelligen
wind_scale = P_wind_prognose_2030 / P_wind




%% Simulate different years
jaren = [2022; 2030];

for jaar = 1:2



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

    %% Scale production for 2030
    if jaren(jaar) == 2030
        P_solar(jaar) = zon_scale .* P_solar;
        P_wind(jaar) = wind_scale .* P_wind;

        PV_sum_producers = zon_scale .* PV_sum_producers;
        Wind_sum_producers = wind_scale .* Wind_sum_producers;
    end

    renewable_producers = PV_sum_producers + Wind_sum_producers;
    residual_load_curve = Consume_curve - renewable_producers;

    residual_fossil_production = residual_load_curve;
    residual_fossil_production(residual_fossil_production<0) = 0;

    % price_electricity = 21.486.*exp(residual_load_curve.*1e-4) ; % v2: y = 21,486e0,0001x, v1: y = 27,775e4E-05x
    price_electricity = 21.486.*exp(residual_load_curve.*1e-4) - (residual_load_curve<0)*21.486; % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production

    price_electricity_raw = price_electricity;
    price_electricity_only_pos = price_electricity(price_electricity>0);

    price_electricity(price_electricity<0) = 0; % set electricity price to zero when residual load is negative = excess energy




    %% Statistics
    % Production volumes
    Prod_annual = sum(Produce_curve)/1000 % [GWh electricity]
    Cons_annual = sum(Consume_curve)/1000 % [GWh electricity]
    Wind_annual = sum(Wind_sum_producers)/1000
    Solar_annual = sum(PV_sum_producers) / 1000
    Prod_wind_perc = Wind_annual / Prod_annual
    Prod_solar_perc = Solar_annual / Prod_annual

    % Electricity Prices
    Price_avg       =   mean(price_electricity)
    Price_max       =   max(price_electricity)
    Price_min       =   min(price_electricity)
    Price_sigma     =   std(price_electricity)
    Price_zero_hours =  length(find(price_electricity==0))
    Price_subzero_hours =  length(find(price_electricity<0))
    Price_only_pos_avg     =   mean(price_electricity_only_pos)








    %% Plot A - generation
    subplot(subplot_size_x,subplot_size_y,  (1 *subplot_size_y)-subplot_size_y   + jaar) % hoogte * y + rij

    plot(time_array,Produce_curve/1000)
    hold on
    xlabel('Time')
    ylabel('Electrical Power [GW]')
    title('Production')

    title(sprintf('Production - Year: %.0f, Consumption: %.1f TWh, %.1f GW Wind, %.1f GWp PV, %.0f prct Wind, %.0f prct PV', jaren(jaar),Cons_annual/1000, P_wind(jaar)/1000  ,P_solar(jaar)/1000,  Prod_wind_perc*100,  Prod_solar_perc*100) )

    % choose time3
    %start_point = 2500; % 7 June 2030
    start_point = 2900; %
    days = 14;
    xlim([time_array(start_point), time_array(start_point+days*24)])
    ylim([0 55])

    % Fossil residual
    area(time_array, (Wind_sum_producers + PV_sum_producers + residual_fossil_production)/1000,'FaceColor','#A2142F') % Purple	'#7E2F8E' , Red ,'#A2142F'

    % Solar
    area(time_array,(Wind_sum_producers + PV_sum_producers)/1000,'FaceColor','#EDB120') % Yellow
    % area(time_array, residual_load_curve+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59})

    % Wind: as last foreground color
    area(time_array,(Wind_sum_producers)/1000,'FaceColor','#77AC30') % Grey


    % ik wil graag overshot ook laten zien met stippelijn erboven over, of negatief?
    % negatief: opladen van batterij

    plot(time_array, Consume_curve/1000,'k')


    legend('Production total','Residual load (mainly fossil backup)','PV solar (household+buildings+central)','Wind energy (inland, coastal, and offshore)','Consumption total')
    grid

    % % sort electricity price
    % duration_curve_price = sort(price_electricity,'descend');
    % area(1:8760,duration_curve_price)
    %
    % %histogram(price_electricity)
    % ylabel('Electricity price €/MWh')
    % xlabel('Occurance [hours per year]')
    % grid
    % % set(gca,'YScale','log')
    %
    %
    % %% Duration curve power
    % duration_curve_power = sort(Consume_curve,'descend')/1000;
    % %duration_curve_power((length(duration_curve_power)+1),1) = 0;
    % area(1:length(duration_curve_power),duration_curve_power)
    % ylim([0 max(duration_curve_power)])
    % xlabel('Hours per year')
    % ylabel('Electrical power [GW]')
    % grid
    % legend('Duration curve consumption')
    % title('Whole year')



    %% Plot B - Residual load
    subplot(subplot_size_x,subplot_size_y,     (2 *subplot_size_y)-subplot_size_y   + jaar) % x*sub_y bepaalt row, sub_y +x bepaalt welke column
    %plot(time_array,Produce_curve)
    %hold on
    %plot(time_array, Consume_curve)
    plot(time_array,residual_load_curve/1e3,'Color','#EDB120') % Why is the residual load sometimes negative? thus more renwable producers than consumers? or renewables that do negative power?
    xlabel('Time')
    ylabel('Electrical Power [GW]')
    grid
    title('Residual load')
    ylim([-30 30])
    xlim([time_array(start_point), time_array(start_point+days*24)])
    % yyaxis right
    % plot(time_array,price_electricity)
    % ylabel('€/MWh')
    %legend('Production total','Consumption total','Residual load','Electricity price')
    legend('Residual load')




    %% Plot C
    subplot(subplot_size_x,subplot_size_y,  (1 *subplot_size_y)-subplot_size_y   + jaar) % hoogte * y + rij

    % plot(time_array,Produce_curve)
    % hold on
    % plot(time_array, Consume_curve)
    % xlabel('Time')
    % ylabel('Electrical Power [MW]')
    % title('Summer week')
    %
    % % choose time
    % start_point = 200; % 7 June 2030
    % days = 300;
    % xlim([time_array(start_point), time_array(start_point+days*24)])
    %
    % % Solar (everything):
    % area(time_array, residual_load_curve+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59}+merit_order_ETM_rawimport{:,44}+merit_order_ETM_rawimport{:,2}+merit_order_ETM_rawimport{:,61}) % wind + stacked: large scale solar
    %
    % % Fossil as second foreground
    % area(time_array, Consume_curve - PV_sum_producers)
    % % area(time_array, residual_load_curve+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59})
    %
    % % Wind: as last foreground color
    % area(time_array, Wind_sum_producers)
    % legend('Production total','Consumption total','PV solar (household+buildings+central)','Residual load (mainly fossil backup)','Wind energy (inland, coastal, and offshore)')
    % grid
    %
    %

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


    %% Plot C - Electricity Price
    subplot(subplot_size_x,subplot_size_y,  ( 3 *subplot_size_y)-subplot_size_y  + jaar) % hoogte +  rij * x

    % plot(time_array,Consume_curve)
    % hold on
    % plot(time_array,residual_load_curve)
    % xlim([time_array(start_point), time_array(start_point+days*24)])

    %ylabel('Electrical Power [MW]')
    %yyaxis right

    plot(time_array,price_electricity,'Color','#D95319')
    hold on
    legend('Day-ahead hourly Electricity price based on merit-order')
    xlabel('Time')
    ylabel('€/MWh')
    grid
    xlim([time_array(start_point), time_array(start_point+days*24)])
    title('Electricity price')
    ylim([0 150])
    %     %
    %     % Electricity Prices
    %     Price_avg       =   mean(price_electricity)
    %     Price_max       =   max(price_electricity)
    %     Price_min       =   min(price_electricity)
    %     Price_sigma     =   std(price_electricity)
    %     Price_zero_hours =  length(find(price_electricity==0))
    %     Price_subzero_hours =  length(find(price_electricity<0))
    %     Price_only_pos_avg     =   mean(price_electricity_only_pos)





    %% Plot D - histogram electricity price
    subplot(subplot_size_x,subplot_size_y,  ( 4 *subplot_size_y)-subplot_size_y  + jaar) % hoogte +  rij * x

    h_elec = histogram(price_electricity);
    h_elec.BinWidth = 5;

    grid
    xlim([0 150])
    xlabel('Electricity price [€/MWh]')
    ylabel('Occurance [hours per year]')
    title('Probability distribution of Electricity price')
    title(sprintf('Electricity prices - avg: %.1f, Std: %.1f, max: %.1f, avg fossil price: %.1f €/MWh',Price_avg, Price_sigma, Price_max, Price_only_pos_avg ) )



end

% Save figure as pdf
% save_fig(h0,'Lipton_PDF_v4_2');     % uses minimized edge borders

% Save figure as png
print -dpng -r300 Lipton_v4_2


%% Find names of largest producers
for j = 1:find_top
    Producer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top(j)}));   % first column: name
    Consumer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_cons_pos_top(j)}));   % first column: name
end
Produced_MWh_year = merit_prod_max_top';
Consumed_MWh_year = merit_cons_max_top';
Largest_prod_cons_table = table(Producer_type_name,Produced_MWh_year,merit_prod_pos_top',Consumer_type_name,Consumed_MWh_year,merit_cons_pos_top');










