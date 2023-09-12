bat1 = [3.6 6.4 10 14.4 18 32 50 72]; % Capacity of battery
cost1 = [400 500 600 700 900 1200 1500 1800]; % Cost of bat1 per respective capacity
NCap = linspace(1,100,50); % The capacities that I want to interpolate and extrapolate the cost to

cost1ext = interp1(bat1, cost1, NCap,'linear', 'extrap');

figure(1)
plot(NCap, cost1ext, '--r')
hold on
plot(bat1, cost1)
hold off
grid
xlabel('Capacity')
ylabel('Cost')