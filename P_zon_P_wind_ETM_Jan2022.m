%% Zon ETM
P_zon_op_huizen = 5610; % [MW]
P_zon_gebouwen = 6248;
P_zonneparken = 5500;
P_zon_op_zee = 0.0;
% Zon in jan 2022: ~10 GWp solar totaal

P_zon_sum_2022 = P_zon_op_huizen + P_zon_gebouwen + P_zonneparken + P_zon_op_zee;

%P_zon_prognose_2030 = 33000; % [MW] laag scenario - als er veel grid congestie is - eprijs dempt flink - curtailment issues in overheidsregeling
P_zon_prognose_2030 = 46200; % [MW] hoog scenario - pv cost down

zon_scale = P_zon_prognose_2030 / P_zon_sum_2022 %[-]

% er zit 15/20/25% piek vermogens verlaging std case in, voor huis, gebouw, park, resp.


%% Wind ETM
P_wind_land = 3300; % [MW] geinstalleerde capaciteit volgens ETM base scenario
P_wind_kust = 1150;
P_wind_zee = 6700;

P_wind_sum_2022 = P_wind_land + P_wind_kust + P_wind_zee;

% Wind in jan 2022: ~10GW wind totaal
P_wind_prognose_2030 = 8800 + 16700; % [MW] laag scenario 8.8GW onshore + 16.7GW offshore
P_wind_prognose_2030 = 8800 + 21300; % [MW] hoog scenario 8.8GW onshore + 21.3GW offshore reeds aangekodigd door overheid, plannen die dit bewerkstelligen

wind_scale = P_wind_prognose_2030 / P_wind_sum_2022 %[-]

