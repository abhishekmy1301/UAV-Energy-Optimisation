% Project workspace cleanup for the transition UAV fidelity example.
% Run during shutdown.
%
% Copyright 2020 -2021 The MathWorks, Inc.

% Clear the figures
delete(findobj('Tag','transitionUAVFidelity'));

% Clear the default parameters
clear FidelityStage;
clear lowFidelityUAV;
clear midFidelityUAV;
clear lowFidHeightResponse;
clear lowFidAirSpeedResponse;
clear lowFidRollResponse;
clear midFidHeightResponse;
clear midFidAirSpeedResponse;
clear midFidRollResponse;


