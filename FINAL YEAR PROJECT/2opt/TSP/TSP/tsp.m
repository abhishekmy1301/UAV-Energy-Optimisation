%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [tour,tourLength] = tsp(position)
% return the optimized tour and its length
% this function shall be used independently from the GUI
% for reference purpose
% position is a list of cities and their coordinates
    tour = [1 2 3];   % initially tour starts with 3 cities
    for i = 4:max(size(position))
        [dist,tour] = tour_insert(tour,position);
    end
    tourLength = dist;
    max1 = 100;
    while(max1)
        n = tour_remove_n(tour);          % generate city to remove
        tournew = tour_remove(tour,n);    % remove city
        [dist,tournew] = tour_insert(tournew,position);
        dist;
        if tourLength>dist  % if distance is lower, update tour, else dont
            tourLength = dist;
            tour = tournew;
        end
        max1 = max1-1;
    end
end