
close all
clear
clc


format short eng
set(groot,'defaultLineLineWidth',2)


%% import methode 1
% onbalans = readtable('Export_01-01-2022_31-01-2022 Balance delta with prices.xls');

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

