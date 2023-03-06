close all
clear
% clc

set(groot,'defaultLineLineWidth',2)

%% Import 12 day subset of data
load("snip_workspace.mat")


tic




%% Old residual load based merit order pricing
% v5 price based on gas and CO2 price
gas = 200; % €/MWh thermal for low and high scenario
plant_eff = 0.55; % [-] thermal to elec eff
gasplant_fuel_cost = gas./plant_eff;
CO2_price = 150; % €/tCO2 emitted in 2022 and in 2030
CO2_emit_gas_plant = 549; % gCO2/kWh = kgCO2/MWh thermal gas
CO2_marginal_cost = CO2_price .* CO2_emit_gas_plant./1000; % €/MWh additional due to CO2 costs
gasplant_marginal_cost = gasplant_fuel_cost + CO2_marginal_cost  % €/MWh electricity delivered marginal costs

fossil_min_price = 59.8; % €/MWh
gasplant_nom_cost_at_residual_load = 20; % GW
bid_incl_exponential = log(gasplant_marginal_cost/fossil_min_price)/(gasplant_nom_cost_at_residual_load*1000);
% exponential price options:
price_electricity = fossil_min_price.*exp(residual_hourly.*bid_incl_exponential); % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production
price_electricity_before_zero = price_electricity;
price_electricity(residual_hourly<0) = 0; % set price to 0 for moments of excess electricity





%% V2G parameters
OBC_power = 11e-3; % [MW] bidirecitonal power transfer capability per EV
n_vehicles = 2.2e6; % [2030]
share_participate_V2G = 0.25; %[-]
n_vehicles_V2G = n_vehicles * share_participate_V2G; % number of vehicles participating in V2G
share_connected_to_charge_pole = 0.20; % 1 out of 5 is connected to charge pole on avg
share_connected_to_charge_pole_V2G = 0.33; % 1 out of 3 of people willing to V2G are connected to charge pole on avg
max_charge_power_all_connected = n_vehicles_V2G * OBC_power; % [MW]
n_vehicles_V2G_connected = n_vehicles_V2G * share_connected_to_charge_pole_V2G; % amount of vehicles that are connected to charging pole AND are willing to do V2G
max_charge_power_inst = n_vehicles_V2G_connected * OBC_power; % [MW]

E_vehicle = 60e-3; %[MWh] storage per vehicle = 65kWh
E_vehicle_V2G_part = 0.50; %[-] --% of SoC is set to be available for V2G
E_vehicle_V2G_fleet = n_vehicles_V2G * E_vehicle * E_vehicle_V2G_part; % [MWh] all V2G vehicles determine max total energy charged since assumed they will connect to charge pole in rolling window and therefore whole battery fleet size can be used - but power is limited by V2G share AND charge pole share
part_V2G_rolling_window = 2.0;
%E_vehicle_V2G_fleet = part_V2G_rolling_window * n_vehicles_V2G * E_vehicle * E_vehicle_V2G_part; % [MWh] all V2G vehicles determine max total energy charged since assumed they will connect to charge pole in rolling window and therefore whole battery fleet size can be used - but power is limited by V2G share AND charge pole share



%% Define cost function = price function
lithium_cost = 119; % €/MWh price in 2025 for pack size
cycle_life = 1850;
deg_cost = lithium_cost./cycle_life.*1000 % degradation cost per MWh throughput of battery in €/MWh



%% Optimizatio: cost function - % include degradation cost in cost function
%fun = @(x) -sum( x.*(56.*exp(bid_incl_exponential.*(residual_hourly-x) ))  -abs(x).*deg_cost/2  ); % 

% improved: deg cost only accounted for during charge, this is during negative 'x', then full deg_cost is accounted
%fun = @(x) -sum( (residual_hourly - x > 0).* x.*(56.*exp(bid_incl_exponential.*(residual_hourly-x) ))  + (residual_hourly - x < 0).*-abs(x).*deg_cost  ); % 

% improved: degradation cost to be seperate sum of x times deg cost. fun negative is profit, fmincon tries to minimize. thus deg_cost should be positive!
% fun notes: -sum( (during discharge       ).* x.*56.*exp(bid_incl_exponential.*( remaining fossil ) )  +sum(abs(x).*deg_cost./2); % 
% fun = @(x) -sum( (residual_hourly - x > 0).* x.*56.*exp(bid_incl_exponential.*(residual_hourly-x)) )  +sum(abs(x).*deg_cost./2); % 

% improved: gains in discharge finance is real electricity price minus deg_cost fixed price
fun = @(x) -sum( (residual_hourly - x > 0).* x.* (56.*exp(bid_incl_exponential.*(residual_hourly-x))-deg_cost) ) ; % 


% fun = @(x) -sum( x.*(56.*exp(bid_incl_exponential.*(residual_hourly-x) )) + x(residual_hourly<x).*-(56.*exp(bid_incl_exponential.*(residual_hourly-x) ))    -abs(x).*deg_cost/2  ); % 
% TODO: mogelijk issue met nog niet ingevulde step function
% TODO: hij moet de juiste residual load tijdstap value pakken
% TODO: hij pakt residual load, maar de elec prijs is 0 als residual < x, dit zit er nog niet in verwerkt.

% TODO: test speed with 8760 hours instead of 337 hours.

% ISSUE: does not seem to find optimal deployment, since not max discharge power is used while prices are higher than an hour later
% ISSUE ANSWER: this is due to degradation cost is integtrated into cost function, therefore minimizes power per hour, but this is not optimal, degradation goes per annual energy!
% TODO: use rolling window for optimization, 

% TODO: normalize x so that value is between 0 and 1, same with price.

% step function; in order to include zero hours price when fossil load = 0; or when residual < x

% starting position
%x0 = zeros(length(time),1);
x0_init = sign(residual_hourly) .* 0.01 .* max_charge_power_inst;

c_discharge = max(residual_hourly) ./ max_charge_power_inst;
c_charge = min(residual_hourly) ./ -max_charge_power_inst;

x0_charge = (residual_hourly + 0.8*c_charge .* max_charge_power_inst < 0 ) .* -0.3 .* max_charge_power_inst;
x0_discharge = (residual_hourly > 0 + 0.8*c_discharge .* max_charge_power_inst) .* 0.3 .* max_charge_power_inst;

% x0 = x0_init + x0_charge + x0_discharge;

x0 = x0_init;


%x0 = residual_hourly/max(residual_hourly) .* 0.1 .* max_charge_power_inst;
%x0 = x0_init;


A1 = tril(ones(length(time)));
A2 = -tril(ones(length(time)));
A = [A1; A2];


b1 = ones(length(time),1).*E_vehicle_V2G_fleet; % maximum charge in battery fleet is this [MWh]
b2 = zeros(length(time),1);                     % minimum charge in battery fleet is 0 [MWh].
b = [b1;b2];
% een inv(A) met A = tril van 8760x8760 duurt 3 seconden.

Aeq = [];
beq = [];

lb = ones(length(time),1).*-max_charge_power_inst;        % mogelijk een lb column array maken met overal dezelfde Pmax waarde
ub = ones(length(time),1).*+max_charge_power_inst;         % ub = [0.5,0.8];

options = optimoptions('fmincon', 'PlotFcn', 'optimplotfval'); % 'Display','iter',
%options = optimoptions('fmincon','MaxFunctionEvaluations', 10000, 'PlotFcn', 'optimplotfval')

options.MaxFunctionEvaluations = 4e5;
options.OptimalityTolerance = 0.0015e6;

%h_fmincon = figure('Name','fmincon'); % ,'pos',[300 100 800 300]

[x,fval]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options); % minimize Xt power in MW * price €/MWh = €
grid


toc 
 % now 11 Aug 2022 it takes 50.9 seconds to run with 337 hours of data, and optimalitytolerance = 0.005e6; and start condition x0 = sign(residual_hourly) .* 0.2 .* max_charge_power_inst; 
%% Plot
h0 = figure('Name','Electricity market NL 2030','pos',[100 200 1400 800]);
h1 = subplot(3,1,1);

plot(time,Wind_hourly ,time,load,time,PV_hourly ,time, residual_hourly)
grid
legend('Wind','Load','PV','Residual')
ylabel('Electrical power [MW]')
% TODO: plot residual after V2G deployment


h2 = subplot(3,1,2);
plot(time,price_electricity)
grid
ylabel('Electricity price €/MWh')
legend('2030')




% determine price after V2G activation
fossil_hourly = residual_hourly - x; % [MW] the first residual (shortage) minus V2G discharge is the remaing shortage = fossil residual load
price_electricity_after = fossil_min_price.*exp(fossil_hourly.*bid_incl_exponential); % [€/MWh] and if residual < 0 than €0/MWh if 0 fossil production or negative residual = excess reneawble energly production
price_electricity_after(fossil_hourly<0) = 0; % set price to 0 for moments of excess electricity



hold on
%plot(time,price_electricity_before_zero)
plot(time,price_electricity_after)
ylim([-250 250])

yyaxis right
plot(time,x)
ylabel('V2G discharge [MW]')
ylim([-max_charge_power_inst +max_charge_power_inst])
legend('Price before V2G','Price after V2G','V2G discharge [MW]')




%% Check 

% revenue
revenue_per_vehicle = -fval/n_vehicles_V2G % € V2G raw revenue per vehicle per total time simulated






% check residual load relative to V2G discharge power
h3 = subplot(3,1,3);
plot(time,residual_hourly)
hold on
plot(time,fossil_hourly)
ylabel('Electrical power [MW]')

plot(time,x)

legend('residual before V2G)','residual after V2G: fossil hourly','V2G discharge power')
grid


%% Save figure

linkaxes([h1 h2 h3],'x')

%save_fig(h0,'FuzeTeaOptim_v2')




% TODO: check energy stored, does it stay within bounds

% probleem is dat ik niet de hoop heb, niet om hulp vraag bij thesis, 



%% check energy stored over time
b0 = figure('Name','check plots','pos',[300 50 1400 600]);
b1 = subplot(2,2,1);

stored_energy = cumsum(x);
plot(time,stored_energy/1000)
hold on
plot( [time(1) time(end)] , [E_vehicle_V2G_fleet/1000 E_vehicle_V2G_fleet/1000])
ylabel('Stored energy in [GWh]')
legend('Stored energy','Max stored energy','Location','Southeast')



%% check x buy and sell prices
b2 = subplot(2,2,2);

% check buy and sell price, and check degradation cost
V2G_hourly_revenue = x .*  price_electricity_after;

% = x.*(56.*exp(bid_incl_exponential.*(residual_hourly-x) )) ; % MW * €/MWh = € [per whole fleet]
deg_cost_hourly = -abs(x).*deg_cost/2;

plot(time,V2G_hourly_revenue./ n_vehicles_V2G)
grid
ylabel('€')
hold on
plot(time,deg_cost_hourly./ n_vehicles_V2G)

legend('revenue','degradation cost')

% cost
income = sum(V2G_hourly_revenue);
cost = sum(deg_cost_hourly);

net_income = income + cost;

income_vehicle = income ./ n_vehicles_V2G
net_income_vehicle = net_income ./ n_vehicles_V2G



%% check max power
% b3 = subplot(2,2,3);
% histogram(x)
% legend('max power')
% xlabel('max instantenous power [MW]')
% grid


%% check x0 initial conditions
b3 = subplot(2,2,3);
plot(time,x0_init,time,x0_charge,time,x0_discharge)
hold on
plot(time,x0,'--')
legend('x0 init','charge','discharge','x0 result')


%% check discharge V2G and price
b4 = subplot(2,2,4);
plot(time,x)
hold on
yyaxis right
plot(time,price_electricity_after)
ylim([-250 250])
grid
legend('V2G discharge','Electricity price')


