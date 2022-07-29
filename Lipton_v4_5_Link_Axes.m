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





%% Simulate different years
jaren = [2022; 2030];



%% Ch.1 Energietransitiemodel input

% import 'NL Power usage (Load)' as an CSV file:
merit_order_ETM_rawimport = readtable('merit_order.911586.csv'); % ETM model - column 2 till 76 is .output - column 77 to 182 is .input powers - all in MW
% converter to matrix for datetime and matrix for doubles
time_array_orig = merit_order_ETM_rawimport{:,1};
time_array = time_array_orig + 365*2 +1;
time_array(:,2) = time_array_orig + 365*10+3;
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
Consume_curve_orig = sum(m_o_consumer{:,2:end},2);

Cons_annual_source_TWh = sum(Consume_curve_orig)/1e6; % ETM source for 2019 is: 124.4 TWh/year which is a lot actually.

E_Load_curves = cat(2,Consume_curve_orig);

Cons_CBS_2019 = 113.4; % [TWh] consumption in NL according to CBS - pre corona source: https://www.cbs.nl/nl-nl/nieuws/2021/09/elektriciteitsproductie-stijgt-in-2020-naar-recordhoogte
Load_CBS_2019 = 118.7; % [TWH] load = consumption + distribution losses, source: https://opendata.cbs.nl/statline/#/CBS/nl/dataset/84575NED/table?ts=1652358103456
Prod_CBS_2019 = 117.6; % [TWh] production of electricity in NL according to CBS

Energy_cons_increase_per_year = 0.03;

Load_expected_demand = [Load_CBS_2019*(1+Energy_cons_increase_per_year*(2022-2019)), Load_CBS_2019*(1+Energy_cons_increase_per_year*(2030-2019))] %[TWh/year divided by TWh/year] source: https://open-pilot.overheid.nl/repository/ronl-f997136c-6917-4bbd-a2f0-5933f3067f67/1/pdf/bijlage-eindrapport-v2g-waarde-en-weg-voorwaarts.pdf
Load_scale = Load_expected_demand / Cons_annual_source_TWh

E_Load_curves  = E_Load_curves .* Load_scale;

Load_jaren = [2019, 2022, 2030]
Load_power_GW_min = [min(Consume_curve_orig)/1000 min(E_Load_curves)/1000] % [GW]
Load_power_GW_mean = [mean(Consume_curve_orig)/1000 mean(E_Load_curves)/1000] % [GW] average hourly load in the three years.
Load_power_GW_max = [max(Consume_curve_orig)/1000 max(E_Load_curves)/1000] % [GW]


if 1 == 2 % histogram van Consumption power NL in 2030 - tussen 12 en 26 GW
    h_hist_load = figure();
    histogram(E_Load_curves(:,1)./1000)
    hold on
    histogram(E_Load_curves(:,2)./1000)
    xlabel('Consumption power in [GW]')
    ylabel('Hourly occurances per year')
    legend(sprintf('Year: 2022, Consumption: %.1f TWh, Mean load: %.1f GW',sum(E_Load_curves(:,1))/1e6, mean(E_Load_curves(:,1)/1000)),sprintf('Year: 2030, Consumption: %.1f TWh, Avg power %.1f GW',sum(E_Load_curves(:,2))/1e6, mean(E_Load_curves(:,2)/1000)) )
    grid
    %title('Histogram of hourly consumption')
    %print -dpng -r300 Histogram_consume_2022_2030
    save_fig(h_hist_load,'Histogram_consume_2022_2030');
end



%% Calculate installed power per type of generator:
% Import production capacities
production_parameters = readtable('production_parameters.911586.csv'); %

%production_parameters.installed_power = production_parameters{:,2} .* production_parameters{:,3};
production_parameters.installed_power = production_parameters{:,"number_of_units"} .* production_parameters{:,"electricity_output_capacity_MW_"};

% Sum up the renewable contributors
P_solar_installed_bron = sum( production_parameters{[14,82,103],"installed_power"} ); % this is shit, row number changes with export ETM, hopefully will not change again, otherwise have to revise script to find based on text string, not on row number
% 17 mei: 11.11 GW PV total van ETM bron

% pv households = 120
% pv buildings = 14
% pv solar parks = 95

P_wind_installed_bron = sum( production_parameters{95:97,"installed_power"} );
% 17 mei: 7.76 GW wind total

% wind onshore 111
% wind coastal 110
% wind offshore 112

% Scale solar and wind to:
%P_zon_prognose_2030 = 33000; % [MW] laag scenario - als er veel grid congestie is - prijs dempt flink - curtailment issues in overheidsregeling
P_zon_2022_April = 14800; %[MW]
P_zon_prognose_2030 = 46200; % [MW] hoog scenario - pv cost down
P_zon_installed_array = [P_zon_2022_April, P_zon_prognose_2030];
zon_scale = P_zon_installed_array / P_solar_installed_bron
% 17 mei: 11.11 GW PV total van ETM bron

P_wind_2022_April = 5300 + 2460; %[MW] 7.76 GW currently onshore + offshore, ratio = 68% wind = onshore (dit zorgt voor kleine fout in modellering, sinds er vanuit ETM bron nu meer offshore dan onshore staat ingesteld)
%P_wind_prognose_2030 = 8800 + 21300; % [MW] hoog scenario 8.8GW onshore + 21.3GW offshore reeds aangekodigd door overheid, plannen die dit bewerkstelligen, maar dit kan niet allemaal nuttig ingevoed worden zonder extra verbruik, dus verwacht: extra H2 electrolysers of extra elec industry
P_wind_prognose_2030 = 8800 + 16700; % [MW] laag scenario 8.8GW onshore + 16.7GW offshore
P_wind_installed_array = [P_wind_2022_April, P_wind_prognose_2030];
wind_scale = P_wind_installed_array / P_wind_installed_bron
% 17 mei: 7.76 GW wind total van ETM bron







%% Step 2: Calculate price
% 2 = PV buildings rooftop solar (industry)
% 44 = PV large scale solar
% 57 = coastal wind energy
% 58 = inland wind energy
% 59 = offshore wind energy
% 61 = PV households
% A) subtract no marginal cost from power usage curve: construct residual load curve
PV_sum_prod_hourly_bron = merit_order_ETM_rawimport{:,2}+merit_order_ETM_rawimport{:,44}+merit_order_ETM_rawimport{:,61} ;
Wind_sum_prod_hourly_bron = merit_order_ETM_rawimport{:,57}+merit_order_ETM_rawimport{:,58}+merit_order_ETM_rawimport{:,59} ;

%% Scale production for all years; 2022 and 2030

PV_sum_prod_hourly = zon_scale .* PV_sum_prod_hourly_bron;
Wind_sum_prod_hourly = wind_scale .* Wind_sum_prod_hourly_bron;

renewable_producers = PV_sum_prod_hourly + Wind_sum_prod_hourly;
residual_load_curves = E_Load_curves - renewable_producers; % [MW] Pure residual load - negative means excess renewable energy produced/available

residual_fossil_production = residual_load_curves;
residual_fossil_production(residual_fossil_production<0) = 0; % [MW] only save positive residual load power hours to only save when fossil is contributing

% v1 price:
% price_electricity = 21.486.*exp(residual_load_curve.*1e-4) ; % v2: y = 21,486e0,0001x, v1: y = 27,775e4E-05x
%price_electricity = 21.486.*exp(residual_load_curves.*1e-4) - (residual_load_curves<0)*21.486; % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production

%v4 price:
% price for 2022 - 2023 non crisis prices used: y = 57,841*exp(0,0516*x
%price_electricity(:,1) = 57.841.*exp(residual_load_curves(:,1).*0.0516e-3) - (residual_load_curves(:,1)<0)*57.841; % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production


% linear merit order price curve according to ETM model, see excel for fit - for 2030:
%price_electricity(:,2) = 7.168/1000.*residual_fossil_production(:,2) + 59.8; % [€/MWh] since using fossil production load only, price is €0 when excess energy - V2G load not taken into account now
%price_electricity(residual_fossil_production==0) = 0; % set price to 0 for moments of excess electricity

% v5 price based on gas and CO2 price
gas = [40; 80]; % €/MWh thermal for low and high scenario
plant_eff = 0.55; % [-] thermal to elec eff
gasplant_fuel_cost = gas./plant_eff;
CO2_price = [80; 150]; % €/tCO2 emitted in 2022 and in 2030
CO2_emit_gas_plant = 549; % gCO2/kWh = kgCO2/MWh thermal gas
CO2_marginal_cost = CO2_price .* CO2_emit_gas_plant./1000; % €/MWh additional due to CO2 costs
gasplant_marginal_cost = gasplant_fuel_cost + CO2_marginal_cost  % €/MWh electricity delivered marginal costs

fossil_min_price = 59.8; % €/MWh
gasplant_nom_cost_at_residual_load = 20; % GW
bid_incl_exponential = log(gasplant_marginal_cost/fossil_min_price)/(gasplant_nom_cost_at_residual_load*1000);
% exponential price options:
price_electricity(:,1) = fossil_min_price.*exp(residual_load_curves(:,1).*bid_incl_exponential(1)); % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production
price_electricity(:,2) = fossil_min_price.*exp(residual_load_curves(:,2).*bid_incl_exponential(2)); % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production
price_electricity(residual_load_curves<0) = 0; % set price to 0 for moments of excess electricity


% make something that if no fossil prod; then elec price is 0.
% if fossil is needed; start price at 59.8

%price_electricity_raw = price_electricity;
% price_electricity_only_pos = price_electricity(price_electricity>0); % has a bug when using array as input, 12918x1 double values instead of expected 8760x2

%price_electricity(price_electricity<0) = 0; % set electricity price to zero when residual load is negative = excess energy

%% Curtailment split between PV and Wind, new: split to ratio of PV/Wind excess energy - OR: to ratio of PV/Wind in ratio to demand.

Curtail_ratio = E_Load_curves ./ renewable_producers .* (residual_load_curves<=0); % if residual load is negative, then calculate curtail ratio
Curtail_ratio(isnan(Curtail_ratio))=0; % a bug occurs when 0 renewable is produced with dividing by zero, therefore this solves it.
Curtail_ratio_topped = Curtail_ratio;
Curtail_ratio_topped = Curtail_ratio_topped + 1 .* (Curtail_ratio==0);

if 1 == 2
    plot(Curtail_ratio)
    hold on
    plot(Curtail_ratio_topped,'--')
    legend('A','B','C','D')
    xlim([0 500])
    ylabel('fraction of renewable energy (Wind+PV) required to meet demand')
end

%% Curtailment of PV
% PV is dominant now, wind is curtailed first, after that PV is curtailed.
%     PV_sum_prod_hourly_curtail = PV_sum_prod_hourly;
%     PV_sum_prod_hourly_curtail(residual_load_curves<0) = Consume_curve(residual_load_curves<0); % limit source to max consumption power of country when source can supply more than country
%     PV_sum_prod_hourly_curtail(PV_sum_prod_hourly<Consume_curve) = PV_sum_prod_hourly(PV_sum_prod_hourly<Consume_curve); % limit source to what is available from source during


% PV is curtailed in ratio of excess energy in balance with wind
PV_sum_prod_hourly_curtail = PV_sum_prod_hourly .* Curtail_ratio_topped;


if 1 == 2 % plot om te controleren:
    h_pv_curt = figure();
    plot(time_array(:,2),PV_sum_prod_hourly(:,2)./1e3,'--')
    hold on
    plot(time_array(:,2),PV_sum_prod_hourly_curtail(:,2)./1e3)
    plot(time_array(:,2),E_Load_curves(:,2)./1e3)
    legend('PV 2030 available','PV 2030 after curtailment','Load curve')
    %legend('zon 2022','zon 2030','zon 2022 curtail','zon 2030 curtail','consume','consume')
    xlim([time_array(2350,2) time_array(2450,2)])
    ylabel('Power [GW]')
    grid
    save_fig(h_pv_curt,'PV_ratio_curtailed');
    print -dpng -r300 PV_ratio_curtailed
end

PV_elec_price = price_electricity; % initiate






%% Curtailment of Wind

%     Wind is dominant:
%     Wind_sum_prod_hourly_curtail = Wind_sum_prod_hourly;
%     Wind_sum_prod_hourly_curtail(residual_load_curves<0) = Consume_curve(residual_load_curves<0); % limit wind to max consumption power of country when wind can supply more than country
%     Wind_sum_prod_hourly_curtail(Wind_sum_prod_hourly<Consume_curve) = Wind_sum_prod_hourly(Wind_sum_prod_hourly<Consume_curve); % limit wind to what is available from wind during

% Wind is curtailed in ratio to excess renewables
Wind_sum_prod_hourly_curtail = Wind_sum_prod_hourly .* Curtail_ratio_topped;


if 1 == 2
    plot(PV_sum_prod_hourly(:,2))
    hold on
    plot(PV_sum_prod_hourly_curtail(:,2))
    plot(E_Load_curves(:,2),'--')
    legend('zon 2030','zon 2030 curtail','consume')
    %legend('zon 2022','zon 2030','zon 2022 curtail','zon 2030 curtail','consume','consume')
    xlim([2560 2660])
    plot(Wind_sum_prod_hourly(:,2))
    plot(Wind_sum_prod_hourly_curtail(:,2),'--')
    plot(Wind_sum_prod_hourly_curtail(:,2)+PV_sum_prod_hourly_curtail(:,2),'.-')
    legend('zon 2030','zon 2030 curtail','consume','wind 2030','wind 2030 curtail','zon+wind curtail')
    grid
    print -dpng -r300 Zon_Wind_ratio_curtailed
end

Wind_elec_price = price_electricity; % initiate


%% PV and Wind electricity prices

PV_elec_price(PV_sum_prod_hourly_curtail == 0) = 0; % set price to zero when PV does not produce, is maybe not necessary since PV production volume is still zero at these instances, thus when multiplying this does not add up, but still handy if non-weighted avg elec price is wanted
PV_revenue_hourly = PV_elec_price .* PV_sum_prod_hourly_curtail; % [€/MWh * MWh] = [€] for every hour
PV_revenue_curt = sum(PV_revenue_hourly); % [€ per year for whole installed base]
PV_avg_revenue_per_MWp = PV_revenue_curt ./ P_zon_installed_array  % P_zon_installed_array
PV_prod_all_annual_volume = sum(PV_sum_prod_hourly); % [MWh]
PV_prod_curt_annual_volume = sum(PV_sum_prod_hourly_curtail); % [MWh]
PV_energy_curtailed_part = (PV_prod_all_annual_volume - PV_prod_curt_annual_volume) ./ PV_prod_all_annual_volume % [ratio]
PV_full_load_factor_avail = PV_prod_all_annual_volume ./ (P_zon_installed_array*365*24) % [hours/year]
PV_full_load_factor_curt = PV_prod_curt_annual_volume ./ (P_zon_installed_array*365*24) % [hours/year]
PV_avg_elec_price_avail = PV_revenue_curt ./ PV_prod_all_annual_volume % [€/MWh]
PV_avg_elec_price_curt = PV_revenue_curt ./ PV_prod_curt_annual_volume % [€/MWh]





%%
if 1 == 2  % plot PV price over year 2022
    plot(time_array,PV_elec_price(:,1))
    ylabel('price [€/MWh]')
    ylim([0 150])
    yyaxis right
    plot(time_array(:,1),PV_sum_prod_hourly_curtail(:,1)./1000)
    ylabel('PV power generated [GW]')
    legend('price PV 2022','power PV 2022')
    %legend('price PV 2022','price PV 2030 curtailed','power PV 2022','power PV 2030')
    ylim([0 30])
    grid
    xlim([time_array(2401) time_array(2545)]) % 11 apri tot 17 april - of: 2 mei t/m 8 mei
    title(sprintf('PV dyn annual volume curtailed: %.1f perc, PV avg worth: %.1f €/MWh', 100.*PV_energy_curtailed_part(1), PV_avg_elec_price_avail(1)) )
    print -dpng -r300 PV_price_2022_week_in_April_title
end

if 1 == 2  % plot PV price over year 2030
    plot(time_array(:,2),PV_elec_price(:,2))
    ylabel('price [€/MWh]')
    yyaxis right
    plot(time_array(:,2),PV_sum_prod_hourly_curtail(:,2)./1000)
    ylabel('PV power generated [GW]')
    legend('price PV','power PV')
    %legend('price PV 2022','price PV 2030 curtailed','power PV 2022','power PV 2030')
    ylim([0 40])
    xlim([time_array(2401,2) time_array(2545,2)]) 
    grid
    print -dpng -r300 PV_price_2030_week_in_April
end

if 1 == 2 % PV dynamic curtailment per year
    h = figure;
    plot(time_array(:,2),PV_sum_prod_hourly_curtail(:,2)./1000)
    hold on
    plot(time_array(:,2),PV_sum_prod_hourly(:,2)./1000,'--')
    ylabel('PV Production [GW]')
    grid
    legend('PV 2030 with dynamic curtailment based on residual load NL','PV 2030 no dynamic curtailment (but with 10-15-20% DC/AC oversizing)')
    xlim([time_array(2401,2) time_array(2545,2)])

    title(sprintf('PV dyn annual volume curtailed: %.1f perc, PV avg worth: %.1f €/MWh', 100.*PV_energy_curtailed_part(2), PV_avg_elec_price_avail(2)) )
    %         yyaxis right
    %         plot(time_array,PV_elec_price(:,2))
    %         ylim([0 250])
    % save_fig(h,'PV curtailment 2030 week in April');
    print -dpng -r300 PV_curtailment_2030_week_in_April
end

Wind_elec_price(Wind_sum_prod_hourly_curtail == 0) = 0;
Wind_revenue_hourly = PV_elec_price .* Wind_sum_prod_hourly_curtail;
Wind_revenue_curt = sum(Wind_revenue_hourly); % [€ per year for whole installed base]
Wind_avg_revenue_per_MW = Wind_revenue_curt ./ P_wind_installed_array' % [€/MW]
Wind_prod_all_annual_volume = sum(Wind_sum_prod_hourly); % [MWh]
Wind_prod_curt_annual_volume = sum(Wind_sum_prod_hourly_curtail); % [MWh]
Wind_energy_curtailed_part = (Wind_prod_all_annual_volume - Wind_prod_curt_annual_volume) ./ Wind_prod_all_annual_volume % [ratio]
Wind_full_load_factor_avail = Wind_prod_all_annual_volume ./ (P_wind_installed_array'*365*24) % [hours/year]
Wind_full_load_factor_curt = Wind_prod_curt_annual_volume ./ (P_wind_installed_array'*365*24) % [hours/year]
Wind_avg_elec_price_avail = Wind_revenue_curt ./ Wind_prod_all_annual_volume % [€/MWh]
Wind_avg_elec_price_curt = Wind_revenue_curt ./ Wind_prod_curt_annual_volume % [€/MWh]



%% Statistics
% Production volumes
Prod_annual = sum(Produce_curve)/1000 % [GWh electricity]
Cons_annual = sum(E_Load_curves)/1000 % [GWh electricity]
Wind_annual = sum(Wind_sum_prod_hourly)/1000
Wind_annual_curtail = sum(Wind_sum_prod_hourly_curtail)/1000
Solar_annual = sum(PV_sum_prod_hourly) / 1000
Solar_annual_curtail = sum(PV_sum_prod_hourly_curtail) / 1000
Prod_wind_perc = Wind_annual_curtail / Prod_annual
Prod_solar_perc = Solar_annual_curtail / Prod_annual

% Electricity Prices
Price_avg       =   mean(price_electricity)
Price_max       =   max(price_electricity)
Price_min       =   min(price_electricity)
Price_sigma     =   std(price_electricity)

Price_zero_hours_moments = price_electricity~=0; % puts a zero when value was zero, else puts a one.
Price_zero_hours_moments = Price_zero_hours_moments - 1; % convert ones to zero, and zeros to minus 1.
Price_zero_hours_moments = Price_zero_hours_moments .* -1; % converts minus 1 to +1
Price_zero_hours = sum(Price_zero_hours_moments) % sums all hours that have zero electricity price per year

%Price_zero_hours =  length(find(price_electricity==0)) % does not work for array values
%Price_subzero_hours =  length(find(price_electricity_raw<0)) % does not work for array values

elec_price_pos_location = price_electricity>0;
Price_only_pos_avg  =  mean(price_electricity.*elec_price_pos_location); %

% deze werkt niet want conditional indexing met matrix A(A>0) output in series array, niet in matrix. Daarom als oplossing: vermenigvuldig met conditional indexing, dan blijft het in matrix format A.*(A>0)
% Price_only_pos_avg(1)  =  mean(price_electricity(price_electricity(:,1)>0,1));
% Price_only_pos_avg(2)  =  mean(price_electricity(price_electricity(:,2)>0,2)) % this shows the price during fossil production hours - non volume weighted




%% 2022 PV - %% Scatter plot to gain insight in PV time-price-production
if 1 == 2
    sl = 1;
    h2020 = figure();
    sz = PV_elec_price(:,sl);
    sz(sz==0) = 1; % set zeros to one since zeros are not accepted by scatter function
    %sz = PV_revenue_hourly(:,sl) ./ 1e4 ; % [€/hour]
    %sz(c==0) = 1;
    c = PV_revenue_hourly(:,sl)/P_zon_installed_array(sl); % [€/hour/MWp installed]
    c(c==0) = 1;
    scatter(time_array(:,1),PV_sum_prod_hourly(:,sl)/1000,c,sz)
    a = colorbar;
    ylabel(a,'Hourly electricity price [€/MWh]','Rotation',270);
    %ylim(a,[0 120])
    %a.Label.String = 'Power (dB)';
    a.Label.Position(1) = 3.2;
    grid
    ylabel('Power produced [GW]')
    legend('2022 14.8 GWp PV, circle size = produced revenue')
    ylim([0 16])
    jaar = 1;
    title(sprintf('Demand: %.1f TWh, %.1f GW Wind, %.1f GWp PV, Wind gen: %.0fprct, PV gen: %.0fprct', Cons_annual(jaar)/1000, P_wind_installed_array(jaar)/1000  ,P_zon_installed_array(jaar)/1000,  Prod_wind_perc(jaar)*100,  Prod_solar_perc(jaar)*100) )
    xlim([time_array(1,1) time_array(end-24,1)])
    print -dpng -r300 Zon_scatterplot_2022_v4
end

%% 2030 PV - %% Scatter plot to gain insight in PV time-price-production
if 1 == 2
    sl = 2;
    h2030 = figure();
    sz = PV_elec_price(:,sl);
    sz(sz==0) = 1; % set zeros to one since zeros are not accepted by scatter function
    %sz = PV_revenue_hourly(:,sl) ./ 1e4 ; % [€/hour]
    %sz(c==0) = 1;
    c = PV_revenue_hourly(:,sl)/P_zon_installed_array(sl); % [€/hour/MWp installed]
    c(c==0) = 1;
    scatter(time_array(:,2),PV_sum_prod_hourly(:,sl)/1000,c*2,sz)
    a = colorbar;
    ylabel(a,'Hourly electricity price [€/MWh]','Rotation',270);
    %a.Label.String = 'Power (dB)';
    a.Label.Position(1) = 3.2;
    grid
    ylabel('Power produced [GW]')
    legend('circle size = produced revenue')
    ylim([0 50])
    xlim([time_array(1,2) time_array(end-24,2)])
    jaar = 2;
    title(sprintf('Demand: %.1f TWh, %.1f GW Wind, %.1f GWp PV, Wind gen: %.0fprct, PV gen: %.0fprct', Cons_annual(jaar)/1000, P_wind_installed_array(jaar)/1000  ,P_zon_installed_array(jaar)/1000,  Prod_wind_perc(jaar)*100,  Prod_solar_perc(jaar)*100) )
    save_fig(h2030,'Zon_scatterplot_2030_v4');
    print -dpng -r300 Zon_scatterplot_2030_v4
end


%% 2022 PV GWp vs Price vs Curtailment per year
jaar = 1;
if 1 == 2
    scatter(PV_sum_prod_hourly(:,jaar)/1000,PV_elec_price(:,jaar))
    xlabel('PV production [GW]')
    ylabel('PV electricity price [€/MWh]')
    grid
%     hold on
%     yyaxis right
%     histogram(PV_sum_prod_hourly(:,jaar)/1000.*PV_elec_price(:,jaar))
%     ylabel('hours per year')
end

%
PV_sum_prod_hourly_GW = PV_sum_prod_hourly./1000;
PV_GW_vs_price = zeros(floor(max(PV_sum_prod_hourly_GW(:,jaar)))+1,1);
PV_GW_vs_revenue = zeros(floor(P_zon_installed_array(:,jaar)./1000),1);
% sum revenue made per GW interval
for a = 1:length(time_array(:,jaar))
     PV_GW_vs_price(floor(PV_sum_prod_hourly_GW(a,jaar))+1) = PV_GW_vs_price(floor(PV_sum_prod_hourly_GW(a,jaar))+1) + PV_elec_price(a,jaar); %[€/MWh]
     PV_GW_vs_revenue(floor(PV_sum_prod_hourly_GW(a,jaar))+1) = PV_GW_vs_revenue(floor(PV_sum_prod_hourly_GW(a,jaar))+1) + PV_elec_price(a,jaar).*PV_sum_prod_hourly(a,jaar)./P_zon_installed_array(:,jaar) ; %[€/MWh * MW instant / MWp] = [€/year/MWp]
% to be finished later - was nice for a twitter post - goal: analyse PV vs GW vs price with no flex in the market.
% sum all things from figure before in vertical direction.
end

PV_GW_vs_revenue_sum = sum(PV_GW_vs_revenue); %[€/year/MWp]
PV_GW_vs_revenue_cumsum = cumsum(PV_GW_vs_revenue); %[€] rising with each additional PV power GW included

PV_GW_vs_revenue_diff = PV_GW_vs_revenue_sum - PV_GW_vs_revenue_cumsum;
PV_revenue_curt_reduction = PV_GW_vs_revenue_diff
%PV_revenue_curt_reduction = circshift(PV_GW_vs_revenue_diff,1);
%PV_revenue_curt_reduction(1) = PV_GW_vs_revenue_sum;
PV_revenue_curt_relative = 1 - PV_revenue_curt_reduction/PV_GW_vs_revenue_sum;

if 1 == 2  % revenue vs curtailment normalized
    figure()
    plot(PV_revenue_curt_relative)
    xlabel('PV [GW] (lower than 14 means curtailment)')
    ylabel('Normalized Annual revenue in euros relative to max income')
    grid
end 

% if 1 == 2
%     figure()
%     bar(PV_GW_vs_price)
%     %xlim([])
%     xlabel('PV power [GW]')
%     ylabel('Cumulative PV elec price')
%     grid
%     legend('jaar030')
% end

if 1 == 2
    figure()
    bar(PV_GW_vs_revenue)
    %xlim([])
    xlabel('PV power [GW]')
    ylabel('Euro per year per MWp per PV power instance')
    grid
    legend('2020')
end



%% 2030 PV GWp vs Price vs Curtailment per year
jaar = 2;
if 1 == 2
    scatter(PV_sum_prod_hourly(:,jaar)/1000,PV_elec_price(:,jaar))
    xlabel('PV production [GW]')
    ylabel('PV electricity price [€/MWh]')
    grid
%     hold on
%     yyaxis right
%     histogram(PV_sum_prod_hourly(:,jaar)/1000.*PV_elec_price(:,jaar))
%     ylabel('hours per year')
end

%
PV_sum_prod_hourly_GW = PV_sum_prod_hourly./1000;
PV_GW_vs_price = zeros(ceil(max(PV_sum_prod_hourly_GW(:,jaar)))+1,1);
PV_GW_vs_revenue = zeros(ceil(P_zon_installed_array(:,jaar)./1000),1);
% sum revenue made per GW interval
for a = 1:length(time_array(:,jaar))
     PV_GW_vs_price(floor(PV_sum_prod_hourly_GW(a,jaar))+1) = PV_GW_vs_price(floor(PV_sum_prod_hourly_GW(a,jaar))+1) + PV_elec_price(a,jaar); %[€/MWh]
     PV_GW_vs_revenue(floor(PV_sum_prod_hourly_GW(a,jaar))+1) = PV_GW_vs_revenue(floor(PV_sum_prod_hourly_GW(a,jaar))+1) + PV_elec_price(a,jaar).*PV_sum_prod_hourly(a,jaar)./P_zon_installed_array(:,jaar) ; %[€/MWh * MW instant / MWp] = [€/year/MWp]
% to be finished later - was nice for a twitter post - goal: analyse PV vs GW vs price with no flex in the market.
% sum all things from figure before in vertical direction.
end

PV_GW_vs_revenue_sum = sum(PV_GW_vs_revenue); %[€/year/MWp]
PV_GW_vs_revenue_cumsum = cumsum(PV_GW_vs_revenue); %[€] rising with each additional PV power GW included

PV_GW_vs_revenue_diff = PV_GW_vs_revenue_sum - PV_GW_vs_revenue_cumsum;
PV_revenue_curt_reduction = PV_GW_vs_revenue_diff;
%PV_revenue_curt_reduction = circshift(PV_GW_vs_revenue_diff,1);
%PV_revenue_curt_reduction(1) = PV_GW_vs_revenue_sum;
PV_revenue_curt_relative = 1 - PV_revenue_curt_reduction/PV_GW_vs_revenue_sum;

if 1 == 2  % revenue vs curtailment normalized
    figure()
    plot(PV_revenue_curt_relative)
    xlabel('PV [GW] (lower than 46.2 GW means curtailment)')
    ylabel('Normalized Annual revenue in euros relative to max income')
    grid
end 


% if 1 == 2
%     figure()
%     bar(PV_GW_vs_price)
%     %xlim([])
%     xlabel('PV power [GW]')
%     ylabel('Cumulative PV elec price')
%     grid
%     legend('jaar030')
% end

if 1 == 2
    figure()
    bar(PV_GW_vs_revenue)
    %xlim([])
    xlabel('PV power [GW]')
    ylabel('Euro per year per MWp per PV power instance')
    grid
    legend('2030')
end







%% 2022 Wind - %% Scatter plot to gain insight in Wind time-price-production
if 1 == 2
    sl = 1;
    h2022_wind = figure();
    sz = Wind_elec_price(:,sl);
    sz(sz==0) = 1; % set zeros to one since zeros are not accepted by scatter function
    %sz = PV_revenue_hourly(:,sl) ./ 1e4 ; % [€/hour]
    %sz(c==0) = 1;
    c = Wind_revenue_hourly(:,sl)/P_wind_installed_array(sl); % [€/hour/MWp installed]
    c(c==0) = 1;
    scatter(time_array(:,1),Wind_sum_prod_hourly(:,sl)/1000,c,sz)
    a = colorbar;
    ylabel(a,'Hourly electricity price [€/MWh]','Rotation',270);
    %ylim(a,[0 120])
    %a.Label.String = 'Power (dB)';
    a.Label.Position(1) = 3.2;
    grid
    ylabel('Power produced [GW]')
    legend('2022 Wind, circle size = produced revenue')
    ylim([0 14])
    jaar = 1;
    title(sprintf('%.0f, Demand: %.1f TWh, %.1f GW Wind, %.1f GWp PV, Wind gen: %.0fprct, PV gen: %.0fprct', jaren(jaar),Cons_annual(jaar)/1000, P_wind_installed_array(jaar)/1000  ,P_zon_installed_array(jaar)/1000,  Prod_wind_perc(jaar)*100,  Prod_solar_perc(jaar)*100) )

    print -dpng -r300 Wind_scatterplot_2022_v2
end

%% 2030 Wind - %% Scatter plot to gain insight in Wind time-price-production
if 1 == 2
    sl = 2;
    h2030_wind = figure();
    sz = Wind_elec_price(:,sl);
    sz(sz==0) = 1; % set zeros to one since zeros are not accepted by scatter function
    %sz = PV_revenue_hourly(:,sl) ./ 1e4 ; % [€/hour]
    %sz(c==0) = 1;
    c = Wind_revenue_hourly(:,sl)/P_wind_installed_array(sl); % [€/hour/MWp installed]
    c(c==0) = 1;
    scatter(time_array(:,2),Wind_sum_prod_hourly(:,sl)/1000,c*2,sz)
    a = colorbar;
    ylabel(a,'Hourly electricity price [€/MWh]','Rotation',270);
    %a.Label.String = 'Power (dB)';
    a.Label.Position(1) = 3.2;
    grid
    ylabel('Power produced [GW]')
    legend('2030 Wind, circle size = produced revenue')
    ylim([0 1.1*max(Wind_sum_prod_hourly(:,sl)/1000)])
    jaar = 2;
    title(sprintf('%.0f, Demand: %.1f TWh, %.1f GW Wind, %.1f GWp PV, Wind gen: %.0fprct, PV gen: %.0fprct', jaren(jaar),Cons_annual(jaar)/1000, P_wind_installed_array(jaar)/1000  ,P_zon_installed_array(jaar)/1000,  Prod_wind_perc(jaar)*100,  Prod_solar_perc(jaar)*100) )

    print -dpng -r300 Wind_scatterplot_2030_v2
end




%% Analyze residual load curve histograms
if 1 == 2
    h_depth = figure();

    sl = 1;
    rx1 = nexttile;
    %subplot(2,1, 1)
    plot(time_array(:,1),residual_load_curves(:,sl)./1000)
    legend('2022')
    grid
    xlabel('Time')
    ylabel('Residual load [GW]')

    sl = 2;
    rx2 = nexttile;
    %subplot(2,1, 2)
    plot(time_array(:,2),residual_load_curves(:,sl)./1000)
    legend('2030')
    grid
    xlabel('Time')
    ylabel('Residual load [GW]')

    linkaxes([rx1 rx2],'x')

    rx1.XLim = [time_array(3000,1), time_array(3000+7*24,1)];
end



%% Analyze volume depth of zero electricity hours
if 1 == 2
    Residual_excess = residual_load_curves;
    Residual_excess(Residual_excess>=0) = NaN; % delete all values that are zero or higher thus fossil required
    Residual_excess_cumsum = cumsum(Residual_excess);
    rs1 = nexttile;
    r2 = histogram(Residual_excess(:,2)./1000,'BinWidth',0.5);
    hold on
    r1 = histogram(Residual_excess(:,1)./1000,'BinWidth',0.5);
    xlabel('Residual load [GW]')
    ylabel('Occurances [hours/year]')
    legend('2030','2022','Location','Northwest')
    title('Histogram of excess available renewable electricity')
    grid
    
    rs2 = nexttile;

    area(r2.BinEdges(1:length(r2.BinEdges)-1),cumsum(r2.Values))
    hold on
    area(r1.BinEdges(1:length(r1.BinEdges)-1),cumsum(r1.Values))

    xlabel('Residual load [GW]')
    ylabel('Available [hours/year]')
    grid
    legend('2030 CDF','2022 CDF','Location','Northwest')
    jaar = 2;
    title(sprintf('Year: %.0f: Consumption: %.1f TWh, %.1f GW Wind, %.1f GWp PV', jaren(jaar),Cons_annual(jaar)/1000, P_wind_installed_array(jaar)/1000  ,P_zon_installed_array(jaar)/1000) )


    linkaxes([rs1 rs2],'x')
    xlim([min(Residual_excess(:,2)./1000) 0])

    print -dpng -r300 Excess_energy_depth_v2

end





%% initiliaze plot
h0 = figure('Name','Electricity market NL 2030','pos',[0 0 2000 1200]); % width and height start and end points
tiledlayout(4,2)




%% Calculate Storage methods
Storage_unlim = zeros(length(time_array),1);
P_V2G_discharge = Storage_unlim;
P_V2G_charge = Storage_unlim;
Storage_V2G = Storage_unlim;

for jaar = 1:2

    %% initialize arrays:

    b = 0;      % b is temporary parameter to store value of last iteration
    %
    %     for a = 1:(length(time_array)-1)
    %         if residual_load_curve(a) < 0 % thus excess energy, than storage activated
    %             %not possible: Storage_unlim(a) = Storage_unlim(a-1) + -residual_load_curve(a);
    %             Storage_unlim(a) = b + -residual_load_curve(a); % b is storage amount of last iteration
    %             b = Storage_unlim(a); % save the value here for the next iteration
    %         end
    %     end

    %% Unlimited storage algorithm:
    % initiliaze first point:

    residual_load_curve = residual_load_curves(:,jaar);

    if residual_load_curve(1) < 0
        Storage_unlim(1) = -residual_load_curve(1);
    end
    % try-out: make it array compatible
    %     Storage_unlim = -residual_load_curve(1,:);
    %     Storage_unlim(Storage_unlim<0) = 0;

    % for loop over other points - by adding on top of previous point - it resets to zero when more days are zero and starts over again.
    for a = 2:(length(time_array)-1)
        if residual_load_curve(a) < 0 % thus excess energy, than storage activated
            Storage_unlim(a) = Storage_unlim(a-1) + -residual_load_curve(a);
        end
    end



    OBC_power = 11e-3; % [MW] bidirecitonal power transfer capability per EV
    n_vehicles = 2.2e6; % [2030]
    share_participate_V2G = 0.40; %[-]
    n_vehicles_V2G = n_vehicles * share_participate_V2G; % number of vehicles participating in V2G
    share_connected_to_charge_pole = 0.33; % 1 out of 5 is connected to charge pole on avg
    max_charge_power_all_connected = n_vehicles_V2G * OBC_power; % [MW]
    n_vehicles_V2G_connected = n_vehicles_V2G * share_connected_to_charge_pole; % amount of vehicles that are connected to charging pole AND are willing to do V2G
    max_charge_power_inst = n_vehicles_V2G_connected * OBC_power; % [MW]

    E_vehicle = 65e-3; %[MWh] storage per vehicle = 65kWh
    E_vehicle_V2G_part = 0.50; %[-] 50% of SoC is set to be available for V2G
    E_vehicle_V2G_fleet = n_vehicles_V2G * E_vehicle * E_vehicle_V2G_part; % [MWh] all V2G vehicles determine max total energy charged - but power is limited by V2G share AND charge pole share


    for a = 1:(length(time_array)-1)

        if residual_load_curve(a) < 0 % thus excess energy, than storage is charged
            Storage_V2G(a+1) = Storage_V2G(a) + -residual_load_curve(a);
            P_V2G_charge(a+1) = -residual_load_curve(a); % If excess energy charge V2G
            if residual_load_curve(a) < -max_charge_power_inst % check charge power limit
                Storage_V2G(a+1) = Storage_V2G(a) + max_charge_power_inst; % limit storage charging to max power and add energy stored
                P_V2G_charge(a+1) = max_charge_power_inst;
            end
            if Storage_V2G(a) > E_vehicle_V2G_fleet % limit Storage capacity (270GWh is 0.9M EV's with 290kWh V2G volume/year)
                Storage_V2G(a) = E_vehicle_V2G_fleet;
                Storage_V2G(a+1) = E_vehicle_V2G_fleet;
                P_V2G_charge(a) = 0;
                P_V2G_charge(a+1) = 0; % if storage is fully charged - no more charging possible thus 0 MW;
            end

        elseif residual_load_curve(a) > 0 % thus shortage of energy - fossil back up required
            Storage_V2G(a+1) = Storage_V2G(a); % make sure storage is kept neutral if not used
            if Storage_V2G(a) > 0 % make sure storage can not deplete more than was charged before
                Storage_V2G(a+1) = Storage_V2G(a) - residual_load_curve(a); % export of stored energy
                P_V2G_discharge(a) = residual_load_curve(a); % Discharge V2G energy only when Storage Energy remaining is still >0.
                if residual_load_curve(a) > max_charge_power_inst
                    Storage_V2G(a+1) = Storage_V2G(a) - max_charge_power_inst;
                    P_V2G_discharge(a) = max_charge_power_inst;
                end
            else
                Storage_V2G(a) = 0; % make sure storage can not go below zero
                Storage_V2G(a+1) = 0;
            end

        end
    end




%% Plot A - generation
rx1 = nexttile;

    %start_point = 3100 = 11 May
    strtPt = 3200; %
    days = 14;
    endPt = strtPt+days*24;

    % V2G discharge
    area(time_array(strtPt:endPt,jaar), (Wind_sum_prod_hourly(strtPt:endPt,jaar) + PV_sum_prod_hourly(strtPt:endPt,jaar) + residual_fossil_production(strtPt:endPt,jaar))/1000,'FaceColor','#7E2F8E') % Purple = #7E2F8E


    %plot(time_array,Produce_curve/1000)

    hold on
    xlabel('Time')
    ylabel('Electrical Power [GW]')
    title('Production')

    title(sprintf('Year: %.0f, Consumption: %.1f TWh, %.1f GW Wind, %.1f GWp PV, Wind gen: %.0f prct, PV gen: %.0f prct', jaren(jaar),Cons_annual(jaar)/1000, P_wind_installed_array(jaar)/1000  ,P_zon_installed_array(jaar)/1000,  Prod_wind_perc(jaar)*100,  Prod_solar_perc(jaar)*100) )

    % choose time3
    xticks(time_array(strtPt:24:endPt,jaar))
    % xlim([time_array(strtPt,(jaar)), time_array(strtPt+days*24,(jaar))])
    ylim([0 55])




    % Fossil residual
    area(time_array(strtPt:endPt,jaar), (Wind_sum_prod_hourly(strtPt:endPt,jaar) + PV_sum_prod_hourly(strtPt:endPt,jaar) + residual_fossil_production(strtPt:endPt,jaar) - P_V2G_discharge(strtPt:endPt))/1000,'FaceColor','#A2142F') %  - Red = #A2142F

    % V2G charge: Wind + Solar
    area(time_array(strtPt:endPt,jaar),(Wind_sum_prod_hourly(strtPt:endPt,jaar) + PV_sum_prod_hourly(strtPt:endPt,jaar))/1000,'FaceColor','#FF0000') % Bright Red = FF0000

    % Solar (Solar-V2G charge)
    area(time_array(strtPt:endPt,jaar),(Wind_sum_prod_hourly(strtPt:endPt,jaar) + PV_sum_prod_hourly(strtPt:endPt,jaar) - P_V2G_charge(strtPt:endPt))/1000,'FaceColor','#EDB120') % Yellow = EDB120
    % area(time_array, residual_load_curve+merit_order_ETM_rawimport{strtPt:endPt,57}+merit_order_ETM_rawimport{strtPt:endPt,58}+merit_order_ETM_rawimport{strtPt:endPt,59})

    %     % V2G Charge
    %     area(time_array,(Wind_sum_producers + P_V2G_charge)/1000,'FaceColor','#FF0000') % Bright Red

    % Wind: as last foreground color
    area(time_array(strtPt:endPt,jaar),(Wind_sum_prod_hourly(strtPt:endPt,jaar) - P_V2G_charge(strtPt:endPt))/1000,'FaceColor','#77AC30') % Grey


    % ik wil graag overshot ook laten zien met stippelijn erboven over, of negatief?
    % negatief: opladen van batterij

    plot(time_array(strtPt:endPt,jaar), E_Load_curves(strtPt:endPt,jaar)/1000,'k')
    plot(time_array(strtPt:endPt,jaar), (E_Load_curves(strtPt:endPt,jaar) + P_V2G_charge(strtPt:endPt))/1000,'--r')

    legend('V2G discharge','Residual load (mainly fossil backup)','V2G charge','PV solar (residential, commercial, central)','Wind energy (inland, coastal, offshore)','Consumption','Consumption (incl flexible)')
    grid






    %% Plot B - Limited Storage only excess energy stored and fed back
    ax2 = nexttile;

    if 1 == 2
        % Histogram of storage size
        figure
        h1 = histogram(Storage/1e3);
        h1.BinWidth = 5;
        xlim([0 600])
        ylabel('Hours per year')
        hold on
        yyaxis right
        ylabel('Annual coverage [%]')
        h2 = cdfplot(Storage/1e3);
        xlabel('GWh storage capacity')
        ylim([0 1.0])
        legend('Histogram','Cumulative distribution function')
    end


    plot(time_array(strtPt:endPt,jaar),Storage_unlim(strtPt:endPt)/1000)
    hold on
    %plot(time_array,Residual_excess_cumsum/1000)
    plot(time_array(strtPt:endPt,jaar),Storage_V2G(strtPt:endPt)/1000)

    xlabel('Time')
    ylabel('Storage [GWh] positive is charging')
    title('Unlimited Storage charging on excess residual load and discharging on shortage')
    grid
    ylim([0 50])
    legend('Unlimited storage charging on excess energy',sprintf('Energy storage in %.0f V2G EVs with a fleet storage of %.f GWh',n_vehicles_V2G,E_vehicle_V2G_fleet/1000))





    %% Plot C - Electricity Price
    ax3 = nexttile;

    % plot(time_array,Consume_curve)
    % hold on
    % plot(time_array,residual_load_curve)
    % xlim([time_array(start_point), time_array(start_point+days*24)])

    %ylabel('Electrical Power [MW]')
    %yyaxis right


    plot(time_array(strtPt:endPt,jaar),price_electricity(strtPt:endPt,jaar),'Color','#D95319')
    legend('Electricity price')
    xlabel('Time')
    ylabel('€/MWh')
    grid
    xlim([time_array(strtPt,jaar), time_array(endPt,jaar)])
    %title('Electricity prices')
    title(sprintf('Electricity prices - Annual avg: %.1f, sigma: %.1f, Free electricity: %.0f hours/year',Price_avg(jaar), Price_sigma(jaar), Price_zero_hours(jaar) ) )

    ylim([0 250])

    %     plot(time_array,P_V2G_charge/1000)
    %     hold on
    %     plot(time_array,P_V2G_discharge/1000)
    %     legend('Charge','Discharge')
    %     xlabel('Time')
    %     ylabel('V2G charge/discharge power GW')
    %     grid
    %     xlim([time_array(start_point), time_array(start_point+days*24)])
    %     title('V2G charging and discharging')
    %     %ylim([0 150])

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
    h_elec = histogram(price_electricity(:,jaar));
    h_elec.BinWidth = 5;

    grid
    xlim([0 250])
    ylim([0 1200])
    xlabel('Electricity price [€/MWh]')
    ylabel('Occurance [hours per year]')
    title('Probability distribution of Electricity price')
    title(sprintf('Electricity prices - PV dyn curtailed: %.1f perc, PV avg worth: %.1f €/MWh, Wind dyn curtialed: %.1f perc, Wind avg worth: %.1f €/MWh', 100.*PV_energy_curtailed_part(jaar), PV_avg_elec_price_avail(jaar), 100.*Wind_energy_curtailed_part(jaar), Wind_avg_elec_price_avail(jaar) ) )





    linkaxes([rx1 ax2 ax3],'x')
    rx1.XLim = [time_array(strtPt,jaar), time_array(endPt,jaar)];

end




% % nog eens:
%     xlim([time_array(start_point), time_array(start_point+days*24)])
%     ylim([0 55])

% Save figure as pdf
% save_fig(h0,'Lipton_PDF_v4_2');     % uses minimized edge borders

%% Save figure as png

% save_fig(h0,'Lipton_v4_5_28_Juli');

% print -dpng -r300 Lipton_v4_5_28_Juli


%% Find names of largest producers
for j = 1:find_top
    Producer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_prod_pos_top(j)}));   % first column: name
    Consumer_type_name(j,:) = convertCharsToStrings((merit_order_ETM_rawimport.Properties.VariableNames{merit_cons_pos_top(j)}));   % first column: name
end
Produced_MWh_year = merit_prod_max_top';
Consumed_MWh_year = merit_cons_max_top';
Largest_prod_cons_table = table(Producer_type_name,Produced_MWh_year,merit_prod_pos_top',Consumer_type_name,Consumed_MWh_year,merit_cons_pos_top');










