close
% Corrected electricity data for each quarter since 2015-2022
electricityData = [2712, 2076, 1907, 2571, 2745, 2065, 1878, 2601, 2623, 1976, 1891, 2571, 2777, 1996, 1786, 2474, 2638, 1979, 1775, 2527, 2655, 2080, 1884, 2664, 2866, 2207, 1803, 2454, 2679, 2035, 1775, 2292];

% Natural gas data for each quarter since 2015-2022
gasData = [11586, 4412, 2335, 7255, 11035, 4591, 1844, 8830, 10424, 4004, 2208, 8736, 12424, 3702, 1878, 8244, 10080, 4425, 1628, 9122, 10166, 4153, 1984, 9196, 11866, 5641, 1895, 8414, 9939, 3770, 1736, 7168];

% Create the heating data by multiplying each value in gasData by 0.8
heatingData = gasData * 0.8;

% Add heating and electricity data
combinedData = heatingData + electricityData;

% Define the quarters and years
quarters = {'Q1', 'Q2', 'Q3', 'Q4'};
years = 2015:2022;

% Create x-axis labels
xLabels = {};
for y = years
    xLabels{end+1} = num2str(y);
    for q = quarters(2:end)
        xLabels{end+1} = '';
    end
end



% Plotting the combined data along with electricity and gas data on the original figure
figure;
hold on;
plot(electricityData, 'b');
plot(gasData, 'r');
plot(heatingData, 'g');
hold off;

xlabel('Year');
ylabel('Energy (thousand tonnes of oil)');
xticks(1:numel(xLabels));
xticklabels(xLabels);
xtickangle(45);
title('Electricity, Natural Gas, and Heating Data for Each Quarter (2015-2022)');
legend('Electricity', 'Natural Gas', 'Heating');




% Plotting the combined data on a new figure
figure;
plot(combinedData, 'b');

xlabel('Year');
ylabel('Required Electricity (thousand tonnes of oil)');
xticks(1:numel(xLabels));
xticklabels(xLabels);
xtickangle(45);
title('Expected Required Electricity for Each Quarter (2015-2022)');
legend('Electricity');


