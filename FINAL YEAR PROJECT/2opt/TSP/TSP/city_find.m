%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [out] = city_find(tour,city)
% find if city is already in the tour
% return 0 if could not find
% return 1 if found
    n = max(size(tour));
    for i = 1:n
        if (tour(i)==city)
            out = 1; break;
        else
            out = 0;
        end
    end
end