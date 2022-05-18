    43'C - 20'C = 23'K bij 800W/m2

    irr = 1000;
    Ppv_raw = (1050Wp * irr/1000)  % LY ES PV: 1050Wp STC: @1000W/m2 @25'C
    dTcell = irr/800*23; % temp rise due to irradiance
    Tcell_abs = Tamb + dTcell;
    P_pv_temp_comp = Ppv_raw * (1 + (Tcell_abs - 25) * -0.39/100)