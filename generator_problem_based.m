clc
clear
close all


load dispatchPrice; % Get poolPrice, which is the revenue per MWh
bar(poolPrice,.5)
xlim([.5,48.5])
xlabel('Price per MWh at each period')
grid



fuelPrice = 3;
totalFuel = 3.95e4;
nPeriods = length(poolPrice); % 48 periods
nGens = 2; % Two generators
gen = [61,152;50,150]; % Generator 1 low = 61 MW, high = 152 MW
fuel = [427,806;325,765]; % Fuel consumption for generator 2 is low = 325, high = 765
startCost = 1e4; % Cost to start a generator after it has been off



efficiency = gen./fuel; % Calculate electricity per unit fuel use
rr = efficiency'; % for plotting
h = bar(rr);
h(1).FaceColor = 'g';
h(2).FaceColor = 'c';
legend(h,'Generator 1','Generator 2','Location','NorthEastOutside')
ax = gca;
ax.XTick = [1,2];
ax.XTickLabel = {'Low','High'};
ylim([.1,.2])
ylabel('Efficiency')




y = optimvar('y',nPeriods,nGens,{'Low','High'},'Type','integer','LowerBound',0,...
    'UpperBound',1);
z = optimvar('z',nPeriods,nGens,'Type','integer','LowerBound',0,...
    'UpperBound',1);



powercons = y(:,:,'Low') + y(:,:,'High') <= 1;





yFuel = zeros(nPeriods,nGens,2);
yFuel(:,1,1) = fuel(1,1); % Fuel use of generator 1 in low setting
yFuel(:,1,2) = fuel(1,2); % Fuel use of generator 1 in high setting
yFuel(:,2,1) = fuel(2,1); % Fuel use of generator 2 in low setting
yFuel(:,2,2) = fuel(2,2); % Fuel use of generator 2 in high setting

fuelUsed = sum(sum(sum(y.*yFuel)));



fuelcons = fuelUsed <= totalFuel;



w = optimexpr(nPeriods,nGens); % Allocate w
idx = 1:(nPeriods-1);
w(idx,:) = y(idx+1,:,'Low') - y(idx,:,'Low') + y(idx+1,:,'High') - y(idx,:,'High');
w(nPeriods,:) = y(1,:,'Low') - y(nPeriods,:,'Low') + y(1,:,'High') - y(nPeriods,:,'High');
switchcons = w - z <= 0;




generatorlevel  = zeros(size(yFuel));
generatorlevel(:,1,1) = gen(1,1); % Fill in the levels
generatorlevel(:,1,2) = gen(1,2);
generatorlevel(:,2,1) = gen(2,1);
generatorlevel(:,2,2) = gen(2,2); 




revenue = optimexpr(size(y));
for ii = 1:nPeriods
    revenue(ii,:,:) = poolPrice(ii)*y(ii,:,:).*generatorlevel(ii,:,:);
end




fuelCost = fuelUsed*fuelPrice;




startingCost = z*startCost;




profit = sum(sum(sum(revenue))) - fuelCost - sum(sum(startingCost));




dispatch = optimproblem('ObjectiveSense','maximize');
dispatch.Objective = profit;
dispatch.Constraints.switchcons = switchcons;
dispatch.Constraints.fuelcons = fuelcons;
dispatch.Constraints.powercons = powercons;



options = optimoptions('intlinprog','Display','final');



[dispatchsol,fval,exitflag,output] = solve(dispatch,'options',options);




%% examine the solution
subplot(3,1,1)
bar(dispatchsol.y(:,1,1)*gen(1,1)+dispatchsol.y(:,1,2)*gen(1,2),.5,'g')
xlim([.5,48.5])
ylabel('MWh')
title('Generator 1 Optimal Schedule','FontWeight','bold')
subplot(3,1,2)
bar(dispatchsol.y(:,2,1)*gen(2,1)+dispatchsol.y(:,2,2)*gen(2,2),.5,'c')
title('Generator 2 Optimal Schedule','FontWeight','bold')
xlim([.5,48.5])
ylabel('MWh')
subplot(3,1,3)
bar(poolPrice,.5)
xlim([.5,48.5])
title('Energy Price','FontWeight','bold')
xlabel('Period')
ylabel('$ / MWh')



starttimes = find(round(dispatchsol.z) == 1); % Use round for noninteger results
[theperiod,thegenerator] = ind2sub(size(dispatchsol.z),starttimes)




