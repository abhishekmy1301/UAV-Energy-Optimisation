function transitionUAVFidelityHighLevelPath()
%TRANSITIONUAVFIDELITYHIGHLEVELPATH simulates and visualzies path following
%characteristics of medium fidelity plant assuming
%transitionUAVFidelityMidLevelSetup has been run before to setup
%parameters.

% Copyright 2020-2021 The MathWorks, Inc.

pathModel='FidelityPathFollowing';

if ~bdIsLoaded(pathModel)
    open_system(pathModel);
end

%Set FidelityStage as 3 for High Fidelity.
FidelityStage=3;
assignin('base','FidelityStage',FidelityStage);

%set maximum roll to 20 degrees.
midFidelityUAV=mediumFidelitySetup();
midFidelityUAV.maxroll=20;
assignin('base','midFidelityUAV',midFidelityUAV);

sim(pathModel);

%Visualize Trajectory
highStates=States;
visualizeSimStates(highStates,FidelityStage);
end



