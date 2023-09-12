%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [dist,newtour] = tour_mindist(tour,city,n)
% return newtour and dist that has min distance when adding city n to the tour
    out = zeros(max(size(city)),max(size(tour))+1);
    for i = 1:max(size(tour))
        tour1 = zeros(1,max(size(tour)+1));
        for l = 1:i
            tour1(l) = tour(l);
        end
        tour1(i+1) = n;
        for l = i+1:max(size(tour))
            tour1(l+1) = tour(l);
        end
        tour1;
        len(i,1) = tour_len(tour1,city);
        out(i,:) = tour1;
    end
    [dist,pos] = min(len);
    newtour = out(pos,:);
end