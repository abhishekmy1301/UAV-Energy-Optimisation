function transitionUAVFidelityMidLevelValidate()
% TRANSITIONUAVFIDELITYSTEPMIDLEVELVALIDATE finds the step response of high
% fidelity plant model and validates the medium fidelity step response.

% Copyright 2020-2021 The MathWorks, Inc.

%Set FidelityStage as 3 for High Fidelity.
stepModel = 'FidelityStepResponse';
assignin('base','FidelityStage',3);

%Gather medium fidelity data.
if(evalin('base','exist(''midFidHeightResponse'')')&&evalin('base','exist(''midFidAirSpeedResponse'')')...
        && evalin('base','exist(''midFidRollResponse'')'))
    midFidHeightResponse=evalin('base',' midFidHeightResponse');
    midFidAirSpeedResponse=evalin('base',' midFidAirSpeedResponse');
    midFidRollResponse=evalin('base',' midFidRollResponse');
else
    %else load stored data.
    load('midFidelityData');
end


%Open the FidelityStepResponse system
modelname = 'FidelityStepResponse';
if ~bdIsLoaded(modelname)
    open_system(modelname);
end

%Find control block in plant model.
controlBlock = get_param('FidelityStepResponse/Step Control Input','Object');

%Get roll response of high fidelity plant.
controlBlock = get_param('FidelityStepResponse/Step Control Input','Object');
controlBlock.StepControl = 'RollAngle Step Control';
sim(stepModel);
highFidRollResponse=highRollAngle;

%Plot the Roll response.
fig_1=figure('Tag','transitionUAVFidelity');
figure(fig_1);
axs=plot(midFidRollResponse.Time,midFidRollResponse.Data(:),'LineWidth',2,'LineStyle',"--");
hold on;
plot(highFidRollResponse.Time,highFidRollResponse.Data(:),'LineWidth',2);
hold off;
axs.Parent.XLim=[75 80];
legend('Mid Fidelity Roll Response','High Fidelity Roll Response');

%Get height response of high fidelity plant.
controlBlock.StepControl = 'Height Step Control';
sim(stepModel);
highFidHeightResponse=highHeight;

%Plot the height response.
fig_2=figure('Tag','transitionUAVFidelity');
figure(fig_2);
axs=plot(midFidHeightResponse.Time,midFidHeightResponse.Data(:),'LineWidth',2,'LineStyle',"--");
hold on;
plot(highFidHeightResponse.Time,highFidHeightResponse.Data(:),'LineWidth',2);
hold off;
legend('Mid Fidelity Height Response','High Fidelity Height Response');
axs.Parent.XLim=[75 80];
%Get AirSpeed response.
controlBlock.StepControl = 'AirSpeed Step Control';
sim(stepModel);

%Plot the Airspeed response.
fig_3=figure('Tag','transitionUAVFidelity');
figure(fig_3);
highFidAirSpeedResponse=highAirSpeed;
axs=plot(midFidAirSpeedResponse.Time,midFidAirSpeedResponse.Data(:),'LineWidth',2,'LineStyle',"--");
hold on;
plot(highFidAirSpeedResponse.Time,highFidAirSpeedResponse.Data(:),'LineWidth',2);
axs.Parent.XLim=[75 80];
hold off;
legend('Mid Fidelity Speed Response','High Fidelity Speed Response');
end
