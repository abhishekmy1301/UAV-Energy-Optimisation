function transitionUAVFidelityMidLevelResponse()
% TRANSITIONUAVFIDELITYMIDLEVELRESPONSE finds the step response of the mid
% fidelity aircraft and plots the results.The low fidelity step response is
% also plotted for comparision.

% Copyright 2020-2021 The MathWorks, Inc.

stepModel = 'FidelityStepResponse';
if ~bdIsLoaded(stepModel)
    open_system(stepModel);
end

assignin('base','FidelityStage',2);
% Get current workspace settings
midFidelityUAV = mediumFidelitySetup;

%set maximum roll in mid fidelity plant
midFidelityUAV.maxroll=35;
assignin('base','midFidelityUAV',midFidelityUAV);

%Gather low fidelity step response data from base workspace if found.
if(evalin('base','exist(''lowFidRollResponse'')')&&evalin('base','exist(''lowFidHeightResponse'')')...
        && evalin('base','exist(''lowFidAirSpeedResponse'')'))
    lowFidRollResponse=evalin('base','lowFidRollResponse');
    lowFidHeightResponse=evalin('base','lowFidHeightResponse');
    lowFidAirSpeedResponse=evalin('base','lowFidAirSpeedResponse');
else
    %else load stored data.
    load('lowFidelityData');
end

%Get roll response for medium fidelity plant
controlBlock = get_param('FidelityStepResponse/Step Control Input','Object');
controlBlock.StepControl = 'RollAngle Step Control';
sim(stepModel);
midFidRollResponse=midRollAngle;

%Plot the roll response.
fig_1=figure('Tag','transitionUAVFidelity');
figure(fig_1);
plot(midFidRollResponse.Time,midFidRollResponse.Data(:),'LineWidth',2);
hold on;
axs=plot(lowFidRollResponse.Time,lowFidRollResponse.Data(:),'LineWidth',2,'LineStyle',"--");
hold off;
legend('Mid Fidelity Roll Response','Low Fidelity Roll Response');
axs.Parent.XLim=[75 80];

%Get height response
controlBlock.StepControl = 'Height Step Control';
sim(stepModel);
midFidHeightResponse=midHeight;

%Plot the height response.
fig_3=figure('Tag','transitionUAVFidelity');
figure(fig_3);
plot(midFidHeightResponse.Time,midFidHeightResponse.Data(:),'LineWidth',2);
hold on;
axs=plot(lowFidHeightResponse.Time,lowFidHeightResponse.Data(:),'LineWidth',2,'LineStyle',"--");
hold off;
legend('Mid Fidelity Height Response','Low Fidelity Height Response');
axs.Parent.XLim=[75 80];

%Get AirSpeed response.
controlBlock.StepControl = 'AirSpeed Step Control';
sim(stepModel);
midFidAirSpeedResponse=midAirSpeed;

%Plot the AirSpeed response.
fig_3=figure('Tag','transitionUAVFidelity');
figure(fig_3);
plot(midFidAirSpeedResponse.Time,midFidAirSpeedResponse.Data(:),'LineWidth',2);
hold on;
axs=plot(lowFidAirSpeedResponse.Time,lowFidAirSpeedResponse.Data(:),'LineWidth',2,'LineStyle',"--");
axs.Parent.XLim=[75 80];
hold off;
legend('Mid Fidelity AirSpeed Response','Low Fidelity AirSpeed Response');

%Assign Mid Fidelity results to workspace for comparision with High
%Fidelity step response.
assignin('base','midFidHeightResponse',midFidHeightResponse);
assignin('base','midFidAirSpeedResponse',midFidAirSpeedResponse);
assignin('base','midFidRollResponse',midFidRollResponse);
end
