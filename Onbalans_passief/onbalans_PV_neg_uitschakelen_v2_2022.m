
close all
clear
clc


format short eng
set(groot,'defaultLineLineWidth',2)


%% import methode 1
onbalans = readtable('Export_01-01-2022_31-01-2022 Balance delta with prices.xls');

%% Plot methode 1
histogram(onbalans.lowest_price_downward)
hold on
histogram(onbalans.highest_price_upward)
grid
legend('afregelen prijs','opregelen')
title('hele maand Jan 2022 - onbalans prijs Tennet')


%% import methode 2
[numbers, strings, raw] = xlsread('Export_01-01-2022_31-01-2022 Balance delta with prices.xls');

%% Plot methode 2
histogram(numbers(:,10))
hold on
histogram(numbers(:,9))
grid
legend('afregelen prijs','opregelen')
title('hele maand Jan 2022 - onbalans prijs Tennet')


%% Calculate
afregelen = numbers(:,10);
opregelen = numbers(:,9);

afregelen_neg = afregelen(afregelen<0) ;
afregelen_pos = afregelen(afregelen>0) ;
afregelen_pos_avg = mean(afregelen_pos)

afregelen_controlled = afregelen;
afregelen_controlled(afregelen_controlled<0) = 0;
afregelen_controlled_clean = rmmissing(afregelen_controlled);
afregelen_controlled_avg = mean(afregelen_controlled_clean)

afregelen_clean = rmmissing(afregelen);
afregelen_avg = mean(afregelen_clean)

opregelen_clean = rmmissing(opregelen);
opregelen_avg = mean(opregelen_clean)


linear_avg_price = ( length(afregelen_clean).*afregelen_avg + length(opregelen_clean).*opregelen_avg ) / (length(afregelen_clean)+length(opregelen_clean))

PV_pos_controlled_avg_price = ( length(afregelen_clean).*afregelen_controlled_avg + length(opregelen_clean).*opregelen_avg ) / (length(afregelen_clean)+length(opregelen_clean))

elec_besparing_per_MWh = PV_pos_controlled_avg_price - linear_avg_price



%% vermenigvuldig met zon PV profiel

% wat is nodig: een PV profiel op 15min, of nog beter op minuut basis
% dat heb ik vanuit uni meting berlijn?


% deze code werkt niet
%[numbers_PV, strings_PV, raw_PV] = xlsread('household_data_1min_singleindex_filtered.csv');

% deze code import wel, maar er zit geen PV power data in
%PV = readtable('household_data_1min_singleindex_filtered.csv');

PV = readtable('NREL_visitor_parking.csv');





%% Plot PV on minute data:
% data points for one year of second data of PV power:
A = 1197509; % 01-01-2017
B = 1435164; % 31-12-2017

PV_installed_power = 524; % [kWp DC]

PV_datetime_2017 = PV.measdatetime(A:B);
PV_power_2017 = PV.ac_power(A:B)./1000; % [kWp AC]

plot(PV_datetime_2017,PV_power_2017)
grid
ylabel('PV output power [kW AC]')
legend('NREL visitor parking 524 kWp DC')

%% annual generation:
PV_annual_generation = sum(PV_power_2017)/60 % [kWh output per year] - divide by 60 since minute recording and want kWh output per year

kWh_per_kWp_NREL = PV_annual_generation / 524 % 1357 kWh/kWp is a pretty sunny location, NL is more like: 980 kWh/kWp, thus let's scale it

kWh_per_kWp_NL = 980;
scale_factor = kWh_per_kWp_NL / kWh_per_kWp_NREL

PV_power_2017_scaled = PV_power_2017 .* scale_factor;

PV_power_2017_scaled_per_kWp = PV_power_2017_scaled ./ PV_installed_power;  % kW AC output for 1 kWp DC installation on minute data scaled to 980 kWh per kWp per year in NL location, roughly

PV_power_2017_scaled_per_MWp = PV_power_2017_scaled_per_kWp ./ 1000;        % MW AC output

%% Tennet whole year imbalance price

imbalance_2022 = readtable('Tennet_imbalance_01012022_31122022.xlsx');

%%
[numbers2, strings2, raw2] = xlsread('Tennet_imbalance_01012022_31122022.xlsx');


%% date time

t1 = datetime([2017 01 01 00 00 00]);
t = t1 + minutes(0:(365*24*60-1));
t = t'; % datetime array in minutes for whole year 2017

%% Plot imbalance whole year

plot(t,-imbalance_2022.Afregelen)
hold on
plot(t,imbalance_2022.Laagste_prijs_afregelen)

plot(t,imbalance_2022.opregelen)
plot(t,imbalance_2022.Hoogste_prijs_opregelen)

legend('afregel vermogen','afregel prijs','opregel vermogen','opregel prijs')
ylabel('MW of €/MWh')
grid

%% Plot PV erbij in

plot(PV_datetime_2017,PV_power_2017)

legend('afregel vermogen','afregel prijs','opregel vermogen','opregel prijs','PV opwek')


%% Retime

tt = timetable(PV_datetime_2017, PV_power_2017_scaled_per_MWp);

% tt3 = retime(tt, 'minutely', 'fillwithmissing');
tt3 = retime(tt, t, 'fillwithmissing');


%% onbalans prijs array construeren
prijs = [imbalance_2022.Laagste_prijs_afregelen, imbalance_2022.Hoogste_prijs_opregelen];

prijs_hoogste = max(prijs,[],2); % [€/MWh] selecteert hoogste prijs voor elk tijdstip tussen af en op regelen

PV_revenue_array = tt3.PV_power_2017_scaled_per_MWp .* prijs_hoogste; % [€/min]
% PV_revenue_array_clean = rmmissing(PV_revenue_array); % this is not useful since it deletes NaN rows.
PV_revenue_array_clean = PV_revenue_array;
PV_revenue_array_clean(isnan(PV_revenue_array)) = 0; % set NaN values to zero, since income/revenue is zero at these time instances

PV_revenue_annual = sum(PV_revenue_array_clean)/60 % [€/year] for a 1 kWp installation

PV_yield = sum(PV_power_2017_scaled_per_MWp)/60*1000 % [kWh per year]


%% if PV inverter does not export when imbalance < 0
clc 

prijs_hoogste_clean = prijs_hoogste;
prijs_hoogste_clean(isnan(prijs_hoogste)) = 0;

prijs_onbalans_only_pos = prijs_hoogste_clean;
prijs_onbalans_only_pos(prijs_hoogste_clean<0) = 0;

% display std revenue
display(PV_revenue_annual)

PV_revenue_array_only_pos = tt3.PV_power_2017_scaled_per_MWp .* prijs_onbalans_only_pos; % [€/min]
PV_revenue_array_only_pos_clean = PV_revenue_array_only_pos;
PV_revenue_array_only_pos_clean(isnan(PV_revenue_array_only_pos)) = 0;
PV_revenue_annual_only_pos_clean = sum(PV_revenue_array_only_pos_clean)/60 

increase_income = PV_revenue_annual_only_pos_clean - PV_revenue_annual



% op een installatie van 1750 kWp DC
% jaaropbrengst PV = ~1750 MWh
% jaaropbrengst waarde PV = 1750 * €102/MWh (SDE bedrag) = 178k
% onbalans opbrengst std PV profiel = €188.0/MWh
% onbalans opbrengst box PV profiel = €200.6/MWh (+12.7/MWh) = 22.1k extra (+7 % !) (op 315k day ahead inkomsten in 2022 bij avg €180/MWh)

% de regel box stuurt PV omvormer vermogen terug naar 0 of 1 % vermogen bij negatieve onbalans prijs
% dat scheelt dus 22k aan onbalanskosten per jaar.

% note: er is wel een derving van minder PV kWh yield hierdoor, maar dat is meegenomen, netto 7% meer financiele opbrengst
% TODO: hoeveel kWh yield lost? dit is een voordeel want dan loopt SDE langer door en minder snel aan 950 kWh/kWp limiet van SDE






%% Plot 4 - PV met onbalans opbrengsten van PV
plot(t,prijs_hoogste)
hold on
plot(PV_datetime_2017,PV_power_2017_scaled_per_MWp.*1e6)
plot(t,PV_revenue_array_clean.*1000) % [€/kWh] om beter duidelijk te maken

grid
legend('onbalans prijs [€/MWh]','PV generation [kW AC]','[€ opbrengst bij 1 MWp PV]')


%%

bin_breedte = 25;

h1 = histogram(PV_power_2017_scaled_per_MWp.*1e6);
h1.Normalization = 'probability';
h1.BinWidth = bin_breedte;

hold on
h2 = histogram(PV_revenue_array.*1000);
h2.Normalization = 'probability';
h2.BinWidth = bin_breedte;

xlabel('kWp production')
ylabel('minute occurances per year')
grid
% ylim([])
legend('PV power of 1MWp DC','PV revenue on imbalance market')




%%


