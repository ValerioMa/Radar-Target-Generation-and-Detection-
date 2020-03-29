%% Implement CFAR
% Implement a Constant False Allarm Rate (CFAR) pick detector on 1D data,
% using lagging cells on the given noise and target scenario.

%% Init 
clear all; close all force; clc; 

%% Params
% Data_points
Ns = 1000;

%% Generate random noise
s=abs(randn(Ns,1));

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s([100 ,200, 300, 700])=[8 9 4 11];


%% Apply CFAR to detect the targets by filtering the noise.

% 1. Define the following:
T = 12;  % Training Cells
G = 4; % Guard Cells 

% Offset : Adding room above noise threshold for desired SNR 
offset=5;

% Vector to hold threshold values 
threshold_cfar = [];

%Vector to hold final signal after thresholding
signal_cfar = [];

% 2. Slide window across the signal length
for i = 1:(Ns-(G+T))     

    % 2. - 5. Determine the noise threshold by measuring it within the training cells
    idxs = 0:(G-1);
    noise_level = sum(s(i:i+T-1));
    threshold = offset * noise_level/T;
    threshold_cfar = [threshold_cfar, {threshold}];
    
    % 6. Measuring the signal within the CUT
    signal = s(i+T+G);
    
    % 8. Filter the signal above the threshold
    if signal<threshold
        signal = 0;
    end
    
    signal_cfar = [signal_cfar, {signal}];
end


% %  signal_cfar is in advance
% % plot the filtered signal
% figure,hold on;
% plot(s);
% plot (cell2mat(signal_cfar),'g--'); 

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')