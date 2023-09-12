%% Solution for The Travelling Salesman Problem:
%%
% A salesman has to visit a number of cities in a single closed tour.
% He always starts and ends the tour in his home city.
% He has to visit each other city on the tour exactly once.
% The solution attempts to minimise the overall travelling distance.
%% For Consultation Services, please contact:
%% Tran Duc Chung, chung.tranduc89@gmail.com || easy-simulink.com
% Enjoy :)

function [out] = distance(a,b)
% calculate distance between a and b
    out = sqrt((a(1)-b(1))^2 + (a(2)-b(2))^2);
end