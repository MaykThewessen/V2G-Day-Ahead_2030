% save workspace part


%start_point = 3100 = 11 May
strtPt = 3200; %
days = 14;
endPt = strtPt+days*24;


jaar = 2;

time = time_array(strtPt:endPt,jaar);
load = E_Load_curves(strtPt:endPt,jaar);
residual_hourly = residual_load_curves(strtPt:endPt,jaar);
Wind_hourly = Wind_sum_prod_hourly(strtPt:endPt,jaar);
PV_hourly = PV_sum_prod_hourly(strtPt:endPt,jaar);


save("snip_workspace.mat","PV_hourly","Wind_hourly","load","residual_hourly","time");


load("snip_workspace.mat")

