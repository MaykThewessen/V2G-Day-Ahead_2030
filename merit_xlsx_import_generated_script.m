%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/mayk/Documents/GitHub/V2G-repository/etm_merit_18_oct_50eu_gas.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 18-Oct-2022 12:21:46

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A3:G27";

% Specify column names and types
opts.VariableNames = ["Technology", "MeritOrderPosition", "MarginalCosts", "InstalledCapacity", "Availability", "FullLoadHoursFuture", "FullLoadHoursPresent"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Technology", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Technology", "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["MeritOrderPosition", "MarginalCosts", "InstalledCapacity", "Availability", "FullLoadHoursFuture", "FullLoadHoursPresent"], "FillValue", 0);

% Import the data
etm_merit_loaded = readtable("/Users/mayk/Documents/GitHub/V2G-repository/etm_merit_18_oct_50eu_gas.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts

