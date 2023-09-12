function transitionUAVFidelityMidLevelPath()
%TRANSITIONUAVFIDELITYMIDLEVELPATH simulates and visualzies path following
%characteristics of medium fidelity plant assuming
%transitionUAVFidelityMidLevelSetup has been run before to setup
%parameters.

% Copyright 2020-2021 The MathWorks, Inc.

%Set FidelityStage as 2 for Medium Fidelity.
FidelityStage=2;
assignin('base','FidelityStage',FidelityStage);

%set maximum roll to 20 degrees.
midFidelityUAV=mediumFidelitySetup();
% midFidelityUAV=evalin('base','midFidelityUAV');
midFidelityUAV.maxroll=20;
assignin('base','midFidelityUAV',midFidelityUAV);

pathModel='FidelityPathFollowing';

if ~bdIsLoaded(pathModel)
    open_system(pathModel);
end

sim(pathModel);

%Visualize Trajectory
mediumStates=States;
visualizeSimStates(mediumStates,FidelityStage);
end

