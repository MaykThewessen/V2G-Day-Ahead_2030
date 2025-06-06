% Histogram of storage size

%Storage = Storage_unlimited;

h0 = figure
h1 = histogram(Storage/1e3);
h1.BinWidth = 5;

ylabel('Hours per year this size of GWh storage is used')
hold on
yyaxis right
ylabel('Annual coverage [%]')

h2 = cdfplot(Storage/1e3);

xlabel('GWh storage capacity')
ylabel('Cumulative % of hours not used per year [-]')
legend('Histogram of storage capacity vs hours per year','Cumulative distribution function vs storage capacity','Location','Southwest')
title('Grid storage market depth and load hours vs capacity - 2030, 30GW Wind, 46GWp PV')


%% save figure
save_fig(h0,'histogram_of_storage_size');
print -dpng -r300 histogram_of_storage_size

%% zoom in
yyaxis left
ylim([0 300])
xlim([0 550])

save_fig(h0,'histogram_of_storage_size_zoom');
print -dpng -r300 histogram_of_storage_size_zoom

%% set log scale
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([5 2000])
ylim([1 500])

save_fig(h0,'histogram_of_storage_size_log_1GW');
print -dpng -r300 histogram_of_storage_size_log_1GW