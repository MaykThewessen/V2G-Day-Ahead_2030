clc
clear
close all


%% residual load

residual = linspace(0,25e3,200);



%% Gas price
gas = [40; 80]; % €/MWh thermal for low and high scenario
plant_eff = 0.55; % [-] thermal to elec eff
gasplant_fuel_cost = gas./plant_eff;
CO2_price = [80; 150]; % €/tCO2 emitted in 2022 and in 2030
CO2_emit_gas_plant = 549; % gCO2/kWh = kgCO2/MWh thermal gas
CO2_marginal_cost = CO2_price .* CO2_emit_gas_plant./1000; % €/MWh additional due to CO2 costs
gasplant_marginal_cost = gasplant_fuel_cost + CO2_marginal_cost  % €/MWh electricity delivered marginal costs





%% fossil price

% v1 price:
% price_electricity = 21.486.*exp(residual_load_curve.*1e-4) ; % v2: y = 21,486e0,0001x, v1: y = 27,775e4E-05x
%price_electricity = 21.486.*exp(residual_load_curves.*1e-4) - (residual_load_curves<0)*21.486; % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production

%v4 price:
% price for 2022 - 2023 non crisis prices used: y = 57,841*exp(0,0516*x
price_electricity(:,1) = 57.841.*exp(residual.*0.0516e-3) - (residual<0)*57.841; % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production


% linear merit order price curve according to ETM model, see excel for fit - for 2030:
price_electricity(:,2) = 7.168/1000.*residual + 59.8; % [€/MWh] since using fossil production load only, price is €0 when excess energy - V2G load not taken into account now
price_electricity(residual==0) = 0; % set price to 0 for moments of excess electricity


%% auto constructed merit price

fossil_min_price = 59.8; % €/MWh
gasplant_nom_cost_at_residual_load = 20; % GW
bid_inclination = (gasplant_marginal_cost - fossil_min_price) / gasplant_nom_cost_at_residual_load % €/MWh increase per GW residual load increase

%linear options:
%price_electricity(:,3) = bid_inclination(1)/1000.*residual + fossil_min_price; % [€/MWh] since using fossil production load only, price is €0 when excess energy - V2G load not taken into account now
%price_electricity(:,4) = bid_inclination(2)/1000.*residual + fossil_min_price; % [€/MWh] since using fossil production load only, price is €0 when excess energy - V2G load not taken into account now

fossil_min_price = 59.8; % €/MWh
gasplant_nom_cost_at_residual_load = 20; % GW
bid_incl_exponential = log(gasplant_marginal_cost/fossil_min_price)/(gasplant_nom_cost_at_residual_load*1000)

% exponential price options:
price_electricity(:,5) = fossil_min_price.*exp(residual.*bid_incl_exponential(1)); % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production
price_electricity(:,6) = fossil_min_price.*exp(residual.*bid_incl_exponential(2)); % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production

price_electricity(residual==0) = 0; % set price to 0 for moments of excess electricity


%% Plot


plot(residual/1e3,price_electricity)

xlabel('Residual load [GW]')
ylabel('Electricity price €/MWh')
grid
legend('2022','2030','2022','2030','2022','2030','Location','northwest')





