%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [tour_data] = tour_dplot(tour,position)
% return tour data for plot
    n = max(size(tour));
    tour = [tour tour(1)];
    pos = zeros(n+1,2); % cities to be plotted
    for i = 1:n+1
        pos(i,1) = position(tour(i),1);
        pos(i,2) = position(tour(i),2);
    end
    tour_data = pos;
end

