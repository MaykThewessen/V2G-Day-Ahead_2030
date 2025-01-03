%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/mayk/Downloads/elaadnl_open_ev_datasets.xlsx
%    Worksheet: open_transactions
%
% Auto-generated by MATLAB on 15-Jul-2022 15:53:55

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 10);

% Specify sheet and range
opts.Sheet = "open_transactions";
opts.DataRange = "A2:J10001";

% Specify column names and types
opts.VariableNames = ["TransactionId", "ChargePoint", "Connector", "UTCTransactionStart", "UTCTransactionStop", "StartCard", "ConnectedTime", "ChargeTime", "TotalEnergy", "MaxPower"];
opts.VariableTypes = ["double", "string", "double", "datetime", "datetime", "string", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["ChargePoint", "StartCard"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ChargePoint", "StartCard"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "UTCTransactionStart", "InputFormat", "");
opts = setvaropts(opts, "UTCTransactionStop", "InputFormat", "");

% Import the data
tbl = readtable("/Users/mayk/Downloads/elaadnl_open_ev_datasets.xlsx", opts, "UseExcel", false);

%% Convert to output type
TransactionId = tbl.TransactionId;
ChargePoint = tbl.ChargePoint;
Connector = tbl.Connector;
UTCTransactionStart = tbl.UTCTransactionStart;
UTCTransactionStop = tbl.UTCTransactionStop;
StartCard = tbl.StartCard;
ConnectedTime = tbl.ConnectedTime;
ChargeTime = tbl.ChargeTime;
TotalEnergy = tbl.TotalEnergy;
MaxPower = tbl.MaxPower;

%% Clear temporary variables
clear opts tbl