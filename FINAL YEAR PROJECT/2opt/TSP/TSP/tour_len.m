%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [len] = tour_len(tour,city)
% calculate tour lenght
% tour = [1 2 3]: tour 1->2->3->1
% city = list of cities
    len=0;
    tour=[tour 1];  % expand tour to cover end point
    n=max(size(tour));
    city1=[city;city(1,:)];
    for i=1:(n-1)
        len=len+distance(city1(tour(i),:),city1(tour(i+1),:));
    end
end