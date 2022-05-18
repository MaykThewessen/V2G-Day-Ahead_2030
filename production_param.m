clear
close all
clc


%% Import production capacities
production_parameters = readtable('production_parameters.904566.csv'); % 

%% Calculate installed power per type of generator:
%production_parameters.installed_power = production_parameters{:,2} .* production_parameters{:,3};
production_parameters.installed_power = production_parameters{:,"number_of_units"} .* production_parameters{:,"electricity_output_capacity_MW_"};


%% Sum up the renewable contributors
P_solar = sum( production_parameters{[14,95,120],"installed_power"} );
% pv households = 120
% pv buildings = 14
% pv solar parks = 95

P_wind = sum( production_parameters{110:112,"installed_power"} );
% wind onshore 111
% wind coastal 110
% wind offshore 112

%% Scale solar and wind to:
%P_zon_prognose_2030 = 33000; % [MW] laag scenario - als er veel grid congestie is - eprijs dempt flink - curtailment issues in overheidsregeling
P_zon_prognose_2030 = 46200; % [MW] hoog scenario - pv cost down
zon_scale = P_zon_prognose_2030 / P_solar

%P_wind_prognose_2030 = 8800 + 16700; % [MW] laag scenario 8.8GW onshore + 16.7GW offshore
P_wind_prognose_2030 = 8800 + 21300; % [MW] hoog scenario 8.8GW onshore + 21.3GW offshore reeds aangekodigd door overheid, plannen die dit bewerkstelligen
wind_scale = P_wind_prognose_2030 / P_wind

%% 

