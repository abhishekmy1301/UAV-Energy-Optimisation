%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [dist,newtour] = tour_insert(tour,city)
% insert a random city in city to a tour
% return new tour that has min distance
    n = max(size(tour));    % no. of cities in current tour
    n0 = max(size(city));   % no. of total cities
    if (n>=n0)
        disp('Cant add more tour.');
    else
        while(1)
            tmp = randi(n0);
            if (0==city_find(tour,tmp))
                [dist,newtour] = tour_mindist(tour,city,tmp);
                break;
            end
        end
    end
end