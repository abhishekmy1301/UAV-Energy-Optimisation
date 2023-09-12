function transitionUAVFidelityMidLevelSetup()
% TRANSITIONFIDELITYMIDSETUP setups aircraft parameters for the mid
% fidelity plant model. 
% Copyright 2020-2021 The MathWorks, Inc.

% load initial values for aircraft.
midFidelityUAV.initialAirspeed=18;
midFidelityUAV.initialAltitude=50;
midFidelityUAV.Tfinal=150;
midFidelityUAV.initialNorth=0;
midFidelityUAV.initialEast=0;
midFidelityUAV.initialFlightPathAngle=0;
midFidelityUAV.initialHeadingAngle=0;

%Load aerodynamic coefficients.
midFidelityUAV.aircraftCL0=0.38;
midFidelityUAV.aircraftDCLda=18.5;
midFidelityUAV.aircraftS=0.9820;
%Aircraft Mass
midFidelityUAV.aircraftM=8.1647;
%Drag polar characteristics.
midFidelityUAV.A1=0.0070;
midFidelityUAV.Apolar=0.0570;
midFidelityUAV.aircraftCD0=0.0220;
midFidelityUAV.aircraftAR=6.9;

midFidelityUAV.MinThk=0.077;
midFidelityUAV.Thk=2.17;
midFidelityUAV.Tfact=0.200;
midFidelityUAV.Kh=1;
midFidelityUAV.KvD=0;

%Choose discrete (0) or continous controller (1)
midFidelityUAV.AirspeedCont=0;
midFidelityUAV.AltCont=0;
%Choose throttle modelling to be advanced (1) or simple model (0).
midFidelityUAV.ThrustAd=1;


%Assign this midFidelityUAV structure to workspace. 
assignin('base','midFidelityUAV',midFidelityUAV);
end
