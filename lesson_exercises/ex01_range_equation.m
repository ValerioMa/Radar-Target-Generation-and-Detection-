%% RANGE EQUATION Calculation
% Using the Radar Range equation it is possile to design the radar 
% transmitter, receiver, and antenna to have the desired power, gain 
% and noise performance to meet the range requirements. 

%Operating frequency (Hz)
fc = 77.0e9;

%Transmitted power (W)
Pt = 3e-3;

%Antenna Gain (linear)
G =  10000;

%Minimum Detectable Power
Ps = 1e-10;

%RCS of a car
RCS = 100;

%Speed of light
c = 3*10^8;

%TODO: Calculate the wavelength
lambda = c/fc;

%TODO : Measure the Maximum Range a Radar can see. 
numerator = Pt*G^2*lambda^2*RCS;
denominator = Ps*(4*pi)^3;
R = (numerator/denominator)^(1/4);