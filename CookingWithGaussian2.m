% Constants for appliances (Watts)
ovenPower = 3000; % Watts
grillPower = 6000; % Watts
kettlePower = 2200; % Watts
toasterPower = 850; % Watts
microwavePower = 800; % Watts

% Time vector for 24-hour period (241 time steps, 6-minute intervals)
timeVec = linspace(0, 24, 241);
numIntervals = length(timeVec);

% Number of households
households = 1000000;

% Initialize power consumption array
powerConsumption = zeros(1, numIntervals);

% Define meal preparation windows and centers (hours)
mealWindows = [7, 9; 12, 14; 18, 20]; % [start, end] for each meal
mealCenters = mean(mealWindows, 2); % Center times for Gaussian distribution

% Appliance durations in minutes
durations = [12, 132, 126]; % [kettle+toaster+microwave, oven, grill]

% Appliance powers for each meal
appliancePowers = [kettlePower + toasterPower + microwavePower; ovenPower; grillPower];

for i = 1:3 % Loop for breakfast, lunch, and dinner
    % Compute the spread (sigma) for the Gaussian distribution
    sigma = (mealWindows(i, 2) - mealWindows(i, 1)) / 6; % Spread (6-sigma rule for ~99% in window)

    % Determine the duration in time steps
    durationMinutes = durations(i);
    durationIntervals = ceil((durationMinutes / 60) * (numIntervals / 24));

    % Get the valid time window in indices
    startHour = mealWindows(i, 1);
    endHour = mealWindows(i, 2);
    windowStart = find(timeVec >= startHour, 1);
    windowEnd = find(timeVec >= endHour, 1);

    % Generate Gaussian-distributed start times using randn
    centerHour = mealCenters(i);
    startTimes = round((centerHour + sigma * randn(households, 1)) * (numIntervals / 24));

    % Clamp start times to the valid window
    startTimes = max(windowStart, min(windowEnd - durationIntervals + 1, startTimes));

    % Build the consumption array for this meal
    for j = 1:households
        startIndex = startTimes(j);
        endIndex = startIndex + durationIntervals - 1;
        powerConsumption(startIndex:endIndex) = powerConsumption(startIndex:endIndex) + appliancePowers(i);
    end
end

% Convert power consumption to MW
powerConsumptionMWcook = powerConsumption / 1e6; % Watts to MW

% Calculate total daily energy in MWh
totalEnergyMWh = sum(powerConsumptionMWcook) * (6 / 60); % MW to MWh

% Print the total daily energy consumption
fprintf('Total daily energy consumed by cooking: %.2f MWh\n', totalEnergyMWh);

% Plotting the 24-hour power consumption for cooking
figure;
plot(timeVec, powerConsumptionMWcook, 'LineWidth', 1.5);
xlabel('Time (Hours)');
ylabel('Power Consumption (MW)');
title('24-Hour Cooking Power Consumption in MW');
grid on;

% Monthly energy consumption
daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
monthlyEnergyMWhcook = totalEnergyMWh * daysInMonth;

% Overlay Gaussian curves for visualization
hold on;
for i = 1:3
    muIdx = find(timeVec >= mealCenters(i), 1);
    sigmaIdx = round((sigma * (numIntervals / 24)));
    gaussian = exp(-0.5 * ((1:numIntervals) - muIdx).^2 / sigmaIdx^2);
    gaussian = gaussian / max(gaussian) * max(powerConsumptionMWcook) / 5; % Normalize for display
    plot(timeVec, gaussian, '--', 'DisplayName', ['Gaussian ' num2str(i)]);
end
legend('Power Demand', 'Gaussian Breakfast', 'Gaussian Lunch', 'Gaussian Dinner');
hold off;
