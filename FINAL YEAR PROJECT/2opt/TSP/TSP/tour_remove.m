%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [out] = tour_remove(tour,city)
% remove a city from a tour
    n = max(size(tour));  
    out = zeros(1,n-1);
    i = 1;
    j = 1;
    while(j<=n)
        if (city~=tour(j))
            out(i) = tour(j);
            i = i+1;
        end
        j = j+1;
    end
end