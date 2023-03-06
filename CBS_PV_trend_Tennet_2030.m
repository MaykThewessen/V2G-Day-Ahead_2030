close all
clear
clc



% format short eng
set(groot,'defaultLineLineWidth',2)

%%
A = [
2010	21
2011	58
2012	138
2013	363
2014	358
2015	519
2016	609
2017	776
2018	1698
2019	2617
2020	3491
2021	3900    
];

PV = table;
PV.year = A(:,1);
PV.added_power = A(:,2);
PV.cum_power = cumsum(PV.added_power);
PV_cum_growth = (PV.cum_power(2:end) - PV.cum_power(1:end-1)) ./ PV.cum_power(1:end-1);


% PV.year = linspace(2010,2030,2030-2010+1)

%% plotting
close all
h0 = figure('Name','CBS trend PV sigmoid','pos',[1200 900 400 250]);
plot(PV.year,PV.cum_power./1000,'-o')
hold on
plot([2021 2030],[PV.cum_power(12)./1000 19.7],'--o')
plot([2021 2030],[PV.cum_power(12)./1000 30.8],'--o')

xlabel('year')
ylabel('Cumulative installed PV capacity [GWp]')
grid
xlim([2010 2030])



% Sigmoid handmatig
x = [PV.year; 2030];
y = [PV.cum_power; 46.2];

year_start = 2021;
year_end = 2030;
years = linspace(year_start,year_end,year_end-year_start+1);
eind_hoogte = 48.5;
schuiven = 2023;
helling = 4/9.5; % in hoeveel jaar van 15 tot 85%?  = 30perc_2y
Y = eind_hoogte ./ (1+exp(-helling.*(years - schuiven)))
plot(years,Y,'--')
plot(2030,46.2,'o','Color','#7E2F8E')
legend('CBS actual installed numbers','Tennet IP2022 scenario: International Ambitie','Tennet IP2022 scenario: Klimaatakkoord','Tennet IP2022 scenario: Nationale Drijfveer','Location','Northwest')

%save_fig(h0,'CBS_PV_trend_sigmoid_tight');

%% Sigmoid function automatic curve fitting toolbox - does not work correctly yet
x = [PV.year; 2030];
y = [PV.cum_power; 46.2];
 
[xData, yData] = prepareCurveData( x, y );

 % Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2010 21];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

