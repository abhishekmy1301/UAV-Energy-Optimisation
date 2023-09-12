%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [out] = tour_remove_n(tour)
% randomly select a number from the tour to remove
    n = max(size(tour));
    m = max(tour);
    found = 0;
    while(1)
        tmp = randi(m);
        for i = 1:n
            if (tmp==tour(i))
                found = 1;
                break;
            end
        end
        if (1==found)
            out = tmp;
            break;
        end
    end
end