% Origin:       14 July 2022
% Author:       Mayk Thewessen
% Department:   Strategy - Research
% Intent:       Electricity market NL analysis for Vehicle-to-Grid in 2030
% Goal:         Charge pole connectedness analysis - Elaad data

close all
clear
clc

% change to test git


%format short eng
set(groot,'defaultLineLineWidth',2)



arrival    = readtable('distribution-of-arrival.xlsx');      % Dataset 2: Distribution of arrival times on weekdays
duration    = readtable('distribution-of-connecti.xlsx');    % Dataset 4: Distribution of connection time per charging event

% Elaad    = readtable('elaadnl_open_ev_datasets.xlsx');    %

%% Arrival time

weight_private = 1/3;
weight_public = 1/3;
weight_workplace = 1/3;
% check:
weight_total = weight_private + weight_public + weight_workplace;

% kwartier waarden:
arrival_avg = (weight_private .* arrival.private + weight_public .* arrival.public + weight_workplace .* arrival.workplace)/100;

for a = 1:25
    b = (a-1) * 4 + 1;
    if a == 25 % fix a bug in a dirty way, since otherwise b would go to #97, while there are only 96 in the dataset
        b = 96-3;
    end
    arrival_hourly_avg(a) = sum(arrival_avg(b:b+3));
end

hours = 0.5:24.5;

bar(hours,arrival_hourly_avg)
grid
xlabel('time')
ylabel('arrival distribution chance per day per hour')


xticks(hours)  % xticks([round(Vehicles_million,1)])
xlim([0 24])
%xtickformat('%.1f M')
%ytickformat('%.2f')



%% Charge duration

duration_cst = 8;

connectedness = arrival_hourly_avg;







%% elaadnl_open_ev_datasets.xlsx

import_Elaad_open % Run .m import file

%Elaad_open = readtable('elaadnl_open_ev_datasets.xlsx');

%% 
figure

h1 = histogram(ChargeTime,'Normalization','probability','BinWidth',1);

hold on

h2 = histogram(ConnectedTime,'Normalization','probability','BinWidth',1);


grid
xlabel('charge duration [hours]')
ylabel('probability of charge duration [-]')

xlim([0 30])

legend('Charge duration','Connected duration')


%% get hour component make histogram

hour_trans_start = hour(UTCTransactionStart);
hour_trans_stop = hour(UTCTransactionStop);

%trans_duration = UTCTransactionStop - UTCTransactionStart;
%hour_trans_duration = round(trans_duration);
hour_trans_duration = round(ConnectedTime);

% delete duration that is 0, replace with 1, since no zeros are accepted
hour_trans_duration(hour_trans_duration<1) = 1;
% add an hour, so hours are from 1 to 24, since this links with matlab numbering starting from 1 onwards
hour_trans_start = hour_trans_start+1;


h3 = histogram(hour_trans_start,'Normalization','probability','BinWidth',1);

hold on

h4 = histogram(hour_trans_stop,'Normalization','probability','BinWidth',1);

legend('start hour','stop hour')
grid


%% calculate start + duration, stack on top of each other

% Let's construct an array that has 24 columns, each represents an hour of the day.
% this array is filled for every hour a transaction/charge pole connection is taking place, consecutively added in a for loop

Connected = zeros(1,24); % first value means from hour 00:00 to 01:00, so middle is 00:30. There are 24 hours per day.


for c = 1:length(hour_trans_start) % loop over every data point = charge session recorder
    e = hour_trans_start(c);
    
    for d = 1:hour_trans_duration(c) % loop over charge duration for every hour
        
        Connected(e) = Connected(e) + 1; % increase value by one, since charge pole is connected
        e = e + 1;
        
        if e > 24 % go round over the night in day, clocks starts at 1 again.
            e = 1;
        end
    end
end

%% Normalization
Connected_sum = sum(Connected);
Connected_norm = Connected ./ Connected_sum;


%% Plotting
h = figure('Name','Charge pole connectedness','pos',[1100 800 450 300]);

bar(0.5:23.5,Connected_norm.*100)
grid
xlabel('hour of the day')
ylabel('distribution over day of connection to charge pole [%]')
legend('n = 10k random public charge events from 2019 open dataset Elaad')
xlim([0 24])

hours = 0:24;
xticks(hours)

save_fig(h,'Charge_pole_connectedness_v3');     % uses minimized edge borders


%% calculate charge duration, instead of charge pole connectedness

% can not be shown, since only data available of duration of charging, not in relation if this is at start middle or end of session.



























