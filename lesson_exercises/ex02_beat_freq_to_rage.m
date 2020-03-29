%% Calculate the range of four targets knowing the beat frequencies
% Given: 
%        the radar maximum range = 300 m 
%        the range resolution = 1m
%        measured beat frequencies [0 MHz, 1.1 MHz, 13 MHz, 24 MHz].
%% 

%% Given params
beat_freq = [0, 1.1, 13, 24]*1e6; % [hz] beat freqency
range_max  = 300;     % [m] radar's max range
c      = 3*10^8;      % [m/s] speed of lights in meters/sec
d_res  = 1;           % [m] range resolution in meters

%% Find the Bsweep of chirp for 1 m resolution ( d_res = c/(2*Bsweep) )
% d_min = c * dt_sweep  --> to have a resolution of d_min we want a wave
% every dt_sweep
% Bsweep = 1/dt_sweep --> from dt extract frequency
Bsweep = c/(2*d_res); % Bsweep calculation

%% Calculate the chirp time (Ts) based on the Radar's Max Range
% To compute the chirp time we need to set a "security" factor that usually
% range betwheen 5 and 6
factor     = 5.5;  % usually between 5 and 6
% Ts = factor*(max travel time)
Ts         = factor*(2*range_max/c); % factor times of the trip time from maximum range


%% Define the frequency shifts 
calculated_range = c*beat_freq*Ts/(2*Bsweep);
% Explanation
% Ts/Bsweep => is the slope of chirp signale x-axis freq, y-axis time
% beat_t = Ts/Bsweep * beat_freq => knowing the beat freq it is possible to extract the
%                                    time
% d = beat_t * c ==> knowing the propagation time and the speed of the
%                    signal it is possible to estimate the propagation
%                    distance (2 way signal i.e. 2xone way signal)
% d_one_way = d/2

%% Display the calculated range
disp(calculated_range);