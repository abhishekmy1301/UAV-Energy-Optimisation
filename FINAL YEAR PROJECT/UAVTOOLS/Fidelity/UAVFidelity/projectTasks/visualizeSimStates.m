function visualizeSimStates(simStates,fidelity)
%VISUALIZESIMSTATES visualizes simulation states logged by pathFollowing model

%   Copyright 2020-2021 The MathWorks, Inc.

x = simStates.North.Data(:);
y = simStates.East.Data(:);
z = -simStates.Height.Data(:);

position = [x y z];

roll = simStates.RollAngle.Data(:);

% In windless condition, body pitch angle is the sum of flight path
% angle and attack angle. For small fixed wing UAV, attack angle is
% usually controlled by autopilot and regulated to be relatively small.
% For visualization purpose, we approximate pitch angle with flight
% path angle.
if (fidelity==3)
    pitch = -simStates.FlightPathAngle.Data(:);
else
    pitch = simStates.FlightPathAngle.Data(:);
end
% In windless, zero side-slip condition, body yaw angle is the same as
% heading angle.
yaw = simStates.HeadingAngle.Data(:);

rotation = quaternion(eul2quat([yaw pitch roll]));

totalPlanes = 50;
plotIndex = 1:floor(size(position,1)/totalPlanes):size(position,1);
plotTransforms(position(plotIndex,:), rotation(plotIndex,:),'MeshFilePath', 'fixedwing.stl','FrameSize',50,'InertialZDirection','Down');

% set the axis ratios
pbaspect([1 1 1]);
daspect([1 1 1]);

% set the grid and light
grid('on');
light('position', [1 0 1]);

end