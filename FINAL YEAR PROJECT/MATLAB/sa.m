clear;clc;close all;
% Define the problem
N = 10; % number of cities
distances = rand(N,N); % randomly generate a distance matrix


% Set initial temperature and other parameters
T0 = 1000; % initial temperature
coolingRate = 0.95; % cooling rate
numIterations = 1000; % number of iterations
numAttempts = 100; % number of attempts to find new solution at each temperature

% Generate initial solution
currentSolution = randperm(N); % randomly generate an initial tour

% Evaluate objective function
currentDistance = calculateDistance(currentSolution, distances); % calculate the total tour distance

% Loop through iterations
for i = 1:numIterations
    % Cooling schedule
    T = T0 * (coolingRate ^ i);
    
    % Try numAttempts times to find a new solution
    for j = 1:numAttempts
        % Generate a new solution
        newSolution = generateNewSolution(currentSolution);
        
        % Evaluate objective function
        newDistance = calculateDistance(newSolution, distances);
        
        % Calculate delta E
        deltaE = newDistance - currentDistance;
        
        % If new solution is better, accept it
        if deltaE < 0
            currentSolution = newSolution;
            currentDistance = newDistance;
        % If new solution is worse, accept it with probability e^(-deltaE/T)
        else
            acceptanceProbability = exp(-deltaE/T);
            if rand() < acceptanceProbability
                currentSolution = newSolution;
                currentDistance = newDistance;
            end
        end
    end
end

% Display final solution and distance
disp('Final tour:');
disp(currentSolution);
disp('Final distance:');
disp(currentDistance);

% Helper functions

function distance = calculateDistance(solution, distances)
    % Calculates the total tour distance for a given solution and distance matrix
    N = length(solution);
    distance = 0;
    for i = 1:N-1
        distance = distance + distances(solution(i), solution(i+1));
    end
    distance = distance + distances(solution(N), solution(1));
end

function newSolution = generateNewSolution(solution)
    % Generates a new solution by swapping two random cities
    N = length(solution);
    i = randi(N);
    j = randi(N);
    while i == j
        j = randi(N);
    end
    newSolution = solution;
    newSolution([i j]) = newSolution([j i]);
end
