function transitionUAVFidelityLowLevelResponse()
% TRANSITIONUAVFIDELITYLOWLEVELRESPONSE finds the step response of the low
% fidelity aircraft.

% Copyright 2020-2021 The MathWorks, Inc.

stepModel = 'FidelityStepResponse';

if ~bdIsLoaded(stepModel)
    open_system(stepModel);
end

%Set Fidelity Stage to 1 to simulate Low Fidelity.
FidelityStage=1;
open_system(stepModel);

%set initial Cruise conditions
lowFidelityUAV.initialAirSpeed=18;
lowFidelityUAV.initialAltitude=50;

%assign parameters to workspace
assignin('base','FidelityStage',FidelityStage);
assignin('base','lowFidelityUAV',lowFidelityUAV);

%Get roll response
controlBlock = get_param('FidelityStepResponse/Step Control Input','Object');
controlBlock.StepControl = 'RollAngle Step Control';
sim(stepModel);
lowFidRollResponse=smRollAngle;
%Get height response
controlBlock.StepControl = 'Height Step Control';
sim(stepModel);
lowFidHeightResponse=smHeight;
%Get AirSpeed response.
controlBlock.StepControl = 'AirSpeed Step Control';
sim(stepModel);
lowFidAirSpeedResponse=smAirSpeed;

%Assign Low Fidelity results to workspace for comparision with Mid
%Fidelity step response.
assignin('base','lowFidHeightResponse',lowFidHeightResponse);
assignin('base','lowFidAirSpeedResponse',lowFidAirSpeedResponse);
assignin('base','lowFidRollResponse',lowFidRollResponse);
end
