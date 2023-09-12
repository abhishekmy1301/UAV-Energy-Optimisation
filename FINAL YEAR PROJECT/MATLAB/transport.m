electricity = [234, 237, 229, 224, 234, 249, 247, 252, 254, 254, 262, 259, 229, 247, 247, 254, 259, 264, 282, 272, 455, 454, 461, 641, 599, 636, 710, 729, 732, 738, 741, 759, 727, 706, 347, 347, 342, 339, 338, 347, 348, 360, 383, 371, 381, 380, 392, 400, 408, 455, 416, 401];
coal_petroleum_coke = [1254, 1186, 1121, 1123, 1048, 1000, 945, 950, 967, 947, 919, 877, 793, 849, 816, 821, 809, 761, 766, 702, 668, 685, 715, 665, 651, 654, 629, 516, 608, 632, 639, 664, 662, 667, 700, 634, 632, 646, 658, 656, 660, 651, 673, 667, 676, 674, 666, 661, 658, 601, 472, 552];

years = 1970:2021; % Years range from 1970 to 2021

% Plot the data
figure;
hold on;
plot(years, electricity, 'b-', 'LineWidth', 2);
plot(years, coal_petroleum_coke, 'r-', 'LineWidth', 2);
hold off;

% Set plot properties
title('Energy Consumption Trends');
xlabel('Years');
ylabel('Energy (thousand tonnes of oil equivalent)');
legend('Electricity', 'Coal, Petroleum, Coke');

% Adjust the axes limits if needed
 xlim([2004, 2022]);
 ylim([300, 750]);

% Add gridlines
grid on;


years = 2004:2021; % Years range from 2004 to 2021

% Data
electricity = [2, 42221;
               2, 42507;
               2, 42513;
               2, 42884;
               2, 41098;
               2, 39635;
               2, 39159;
               2, 38646;
               2, 38508;
               3, 38177;
               6, 38713;
               8, 39510;
               11, 40429;
               16, 40522;
               21, 39959;
               23, 39146;
               44, 31792;
               74, 35229];

% Filter the data for the years starting from 2004
start_index = find(years == 2004);
electricity_filtered = electricity(start_index:end, :);

% Plot the filtered data
figure;
hold on;
plot(years(start_index:end), electricity_filtered(:, 1), 'b-', 'LineWidth', 2);
plot(years(start_index:end), electricity_filtered(:, 2), 'r-', 'LineWidth', 2);
hold off;

% Set plot properties
title('Energy Consumption Trends (from 2004)');
xlabel('Years');
ylabel('Energy (thousand tonnes of oil equivalent)');
legend('Electricity', 'Petroleum, Bioenergy');

% Add gridlines
grid on;
