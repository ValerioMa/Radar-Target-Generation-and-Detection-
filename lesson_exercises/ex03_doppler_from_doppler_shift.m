%% Compute speed from doppler frequency shifts 
% Calculate thevelocity in m/s of four targets with following doppler 
% frequency shifts: .
% You can use the following parameter values:
%     doppler_shift [3 KHz, 4.5 KHz, 11 KHz, -3 KHz]
%     The radar's operating frequency = 77 GHz
%     The speed of light c = 3*10^8

% The beat frequency is not only related to the range of the target but
% also to its relative radial velocity with respect to the radar 

% The doppler shift is a shift in the received signal frequency due to the
% doppler effecto of the target velocity.

%% Parameters:
c = 3*10^8;         %speed of light
frequency = 77e9;   %frequency in Hz
doppler_shift = [3, 4.5, 11, -3]*1e3; % frequency shift in Hz

%% Calculate the wavelength
wavelength = c/frequency;

%% Calculate the velocity of the targets  fd = 2*vr/lambda
Vr = doppler_shift*wavelength/2;

% phase = dist/lambd
% delta_phase = 2*delta_dist/lambda
% delta_phase = 2*delta_dist*freq/speed

% the delta_phase can be translated in freq 
% delta_freq = delta_phase/delta_time

% Therefore
% delta_freq = delta_phase /delta_time = 2*delta_dist/(lambda*delta_time)
% delta_freq = 2*delta_dist/(lambda*delta_time)
% delta_freq = 2*speed/lambda

%% Display results
disp(Vr);