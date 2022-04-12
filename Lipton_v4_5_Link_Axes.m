% Origin:       Tuesday 8 Feb 2022
% Author:       Mayk Thewessen
% Department:   Strategy - Research
% Intent:       Electricity market NL analysis for Vehicle-to-Grid in 2030


close all
clear
clc

% change to test git


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
tiledlayout(4,2)


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

%% Calculate installed power per type of generator:
% Import production capacities
production_parameters = readtable('production_parameters.904566.csv'); %

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




    %% Calculate Storage methods
    
    % initialize arrays:
    Storage = zeros(length(time_array),1);
    P_V2G_discharge = Storage;
    P_V2G_charge = Storage;
    Storage3 = Storage;

    for a = 1:(length(time_array)-1)
        if residual_load_curve(a) < 0 % thus excess energy, than storage activated
            Storage(a+1) = Storage(a) + -residual_load_curve(a);
        end
    end

    Residual_excess = -residual_load_curve;
    Residual_excess(Residual_excess<0) = 0;
    Residual_excess_cumsum = cumsum(Residual_excess);

    OBC_power = 11e-3; % [MW] bidirecitonal power transfer capability per EV
    n_vehicles = 2.2e6; % [2030]
    share_participate_V2G = 0.4; %[-]
    n_vehicles_V2G = n_vehicles * share_participate_V2G; % number of vehicles participating in V2G
    share_connected_to_charge_pole = 0.2; % 1 out of 5 is connected to charge pole on avg
    max_charge_power_all_connected = n_vehicles_V2G * OBC_power; % [MW]
    max_charge_power_inst = n_vehicles_V2G * OBC_power * share_connected_to_charge_pole; % [MW]

    E_vehicle = 65e-3; %[MWh] storage per vehicle
    E_vehicle_V2G_part = 0.5; %[-]
    E_vehicle_V2G_fleet = n_vehicles_V2G * E_vehicle * E_vehicle_V2G_part; % [MWh]
    

    for a = 1:(length(time_array)-1)

        if residual_load_curve(a) < 0 % thus excess energy, than storage is charged
            Storage3(a+1) = Storage3(a) + -residual_load_curve(a);
            P_V2G_charge(a+1) = -residual_load_curve(a); % If excess energy charge V2G
            if residual_load_curve(a) < -max_charge_power_inst % check charge power limit
                Storage3(a+1) = Storage3(a) + max_charge_power_inst; % limit storage charging to max power and add energy stored
                P_V2G_charge(a+1) = max_charge_power_inst;
            end
            if Storage3(a) > E_vehicle_V2G_fleet % limit Storage capacity (270GWh is 0.9M EV's with 290kWh V2G volume/year)
                Storage3(a) = E_vehicle_V2G_fleet;
                Storage3(a+1) = E_vehicle_V2G_fleet;
                P_V2G_charge(a) = 0;
                P_V2G_charge(a+1) = 0; % if storage is fully charged - no more charging possible thus 0 MW;
            end

        elseif residual_load_curve(a) > 0 % thus shortage of energy - fossil back up required
            Storage3(a+1) = Storage3(a); % make sure storage is kept neutral if not used
            if Storage3(a) > 0 % make sure storage can not deplete more than was charged before
                Storage3(a+1) = Storage3(a) - residual_load_curve(a); % export of stored energy
                P_V2G_discharge(a) = residual_load_curve(a); % Discharge V2G energy only when Storage Energy remaining is still >0.                
                if residual_load_curve(a) > max_charge_power_inst
                    Storage3(a+1) = Storage3(a) - max_charge_power_inst;
                    P_V2G_discharge(a) = max_charge_power_inst;
                end
            else
                Storage3(a) = 0; % make sure storage can not go below zero
                Storage3(a+1) = 0;
            end
            
        end
    end






    %% Plot A - generation
    ax1 = nexttile;

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



    % V2G discharge    
    area(time_array, (Wind_sum_producers + PV_sum_producers + residual_fossil_production)/1000,'FaceColor','#7E2F8E') % Purple = #7E2F8E 

    % Fossil residual
    area(time_array, (Wind_sum_producers + PV_sum_producers + residual_fossil_production - P_V2G_discharge)/1000,'FaceColor','#A2142F') %  - Red = #A2142F 

    % V2G charge: Wind + Solar
    area(time_array,(Wind_sum_producers + PV_sum_producers)/1000,'FaceColor','#FF0000') % Bright Red = FF0000

    % Solar (Solar-V2G charge)
    area(time_array,(Wind_sum_producers + PV_sum_producers - P_V2G_charge)/1000,'FaceColor','#EDB120') % Yellow = EDB120
    % area(time_array, residual_load_curve+merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59})

%     % V2G Charge
%     area(time_array,(Wind_sum_producers + P_V2G_charge)/1000,'FaceColor','#FF0000') % Bright Red

    % Wind: as last foreground color
    area(time_array,(Wind_sum_producers - P_V2G_charge)/1000,'FaceColor','#77AC30') % Grey


    % ik wil graag overshot ook laten zien met stippelijn erboven over, of negatief?
    % negatief: opladen van batterij

    plot(time_array, Consume_curve/1000,'k')
    plot(time_array, (Consume_curve + P_V2G_charge)/1000,'--r')

    legend('--','V2G discharge','Residual load (mainly fossil backup)','V2G charge','PV solar (household+buildings+central)','Wind energy (inland, coastal, and offshore)','Consumption (inflexible)','Consumption (incl flexible)')
    grid
    
    




    %% Plot B - Limited Storage only excess energy stored and fed back
    ax2 = nexttile;
    
    %Histogram of storage size
%         figure
%         h1 = histogram(Storage/1e3);
%         h1.BinWidth = 5;
%         xlim([0 600])
%         ylabel('Hours per year')
%         hold on
%         yyaxis right
%         ylabel('Annual coverage [%]')
%         h2 = cdfplot(Storage/1e3);
%         xlabel('GWh storage capacity')
%         ylim([0 1.0])
%         legend('Histogram','Cumulative distribution function')

    plot(time_array,Storage/1000)
    hold on
    %plot(time_array,Residual_excess_cumsum/1000)
    plot(time_array,Storage3/1000)

    xlabel('Time')
    ylabel('Storage [GWh] positive is charging')
    title('Unlimited Storage charging on excess residual load and discharging on shortage')
    grid
    ylim([-1 200])
    legend('Unlimited storage charging on excess energy',sprintf('Energy storage in %.0f V2G EVs with a fleet storage of %.f GWh',n_vehicles_V2G,E_vehicle_V2G_fleet/1000))





    %% Plot C - Electricity Price
    ax3 = nexttile;

    % plot(time_array,Consume_curve)
    % hold on
    % plot(time_array,residual_load_curve)
    % xlim([time_array(start_point), time_array(start_point+days*24)])

    %ylabel('Electrical Power [MW]')
    %yyaxis right

    plot(time_array,P_V2G_charge/1000)
    hold on
    plot(time_array,P_V2G_discharge/1000)
    legend('Charge','Discharge')
    xlabel('Time')
    ylabel('V2G charge/discharge power GW')
    grid
    xlim([time_array(start_point), time_array(start_point+days*24)])
    title('V2G charging and discharging')
    %ylim([0 150])

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
    ax4 = nexttile;
    h_elec = histogram(price_electricity);
    h_elec.BinWidth = 5;

    grid
    xlim([0 150])
    xlabel('Electricity price [€/MWh]')
    ylabel('Occurance [hours per year]')
    title('Probability distribution of Electricity price')
    title(sprintf('Electricity prices - avg: %.1f, Std: %.1f, max: %.1f, avg fossil price: %.1f €/MWh',Price_avg, Price_sigma, Price_max, Price_only_pos_avg ) )





end

linkaxes([ax1 ax2 ax3],'x')
ax1.XLim = [time_array(start_point), time_array(start_point+days*24)];

% % nog eens:
%     xlim([time_array(start_point), time_array(start_point+days*24)])
%     ylim([0 55])

% Save figure as pdf
% save_fig(h0,'Lipton_PDF_v4_2');     % uses minimized edge borders

% Save figure as png
print -dpng -r300 Lipton_v4_5_fix_wind_overlaps_V2G


%% Find names of largest producers
for j = 1:find_top
    Producer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top(j)}));   % first column: name
    Consumer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_cons_pos_top(j)}));   % first column: name
end
Produced_MWh_year = merit_prod_max_top';
Consumed_MWh_year = merit_cons_max_top';
Largest_prod_cons_table = table(Producer_type_name,Produced_MWh_year,merit_prod_pos_top',Consumer_type_name,Consumed_MWh_year,merit_cons_pos_top');










