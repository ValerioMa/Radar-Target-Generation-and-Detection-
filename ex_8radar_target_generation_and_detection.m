clear all
clc;
% close all force;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%speed of light = 3e8
%%%%%%%%%%%%%%%%%%%%%%%%%%%

speed_of_light   = 3e8;
max_velocity     = 100;
range_resolution = 1;
max_range        = 200;

saw_thoot = true;
if saw_thoot 
    crop_idxs = 150;
else
    crop_idxs = 0;
end
%% User Defined Range and Velocity of target
x = [60; 30];  % state of the target [range, speed]

x0 = x;
x0_str = ['(s0 = ', num2str(x0(1)), ', v0 = ', num2str(x0(2)),')'];

%% FMCW Waveform Generation
% Bandwidth
Bsweep = speed_of_light/(2 * range_resolution);

% Chirp Time
chirp_factor = 5.5;
Tchirp = chirp_factor * 2 * max_range/speed_of_light;

% Slope
slope = Bsweep/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
wavelength = speed_of_light/fc;

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024 + 2*crop_idxs;                  %for length of time OR # of range cells
% In the saw tooth profile the first and last part must be cutted away

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Txf=zeros(1,length(t)); %transmitted freq
Rxf=zeros(1,length(t)); %received freq
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 
dt = t(2) - t(1); % we have a constant dt
H = [1, dt ; 0, 1];
for i=1:length(t)         
    
    
    %For each time stamp update the Range of the Target for constant velocity. 
    x = H*x;
    r_t(i) = x(1);
    td(i) = 2*x(1)/speed_of_light; % approximate delay (we consider the car stationary)
    deltaT = t(i) - td(i);
    

    %For each time sample we need update the transmitted and received signal.
    if saw_thoot % depends on the selected profile
        [Tx(i), Txf(i)] = sawtoothWave(t(i), fc, Tchirp, Bsweep);
        [Rx(i), Rxf(i)] = sawtoothWave(deltaT, fc, Tchirp, Bsweep);
    else
        Txf(i) = fc + slope * t(i)/2; 
        Rxf(i) = fc + slope * deltaT/2; 
        Mix(i) = cos( 2 * pi * (Txf(i) * t(i)));    
        Rx(i) = cos( 2 * pi * (Rxf(i) * deltaT));
    end
    

    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) * Rx(i);
    %cos(2*pi*(2*Slope*x(1)/c*t(i) + 2*fc*x(2)/c*t(i)))
    
end

%% Reshape the vector into Nr*Nd array.
% Nr and Nd here would also define the size of Range and Doppler FFT respectively.
Txf_t = reshape(Txf, Nr, Nd);
Rxf_t = reshape(Rxf, Nr, Nd);
Tx_t  = reshape(Tx, Nr, Nd);
Rx_t  = reshape(Rx, Nr, Nd);
Mix_t = reshape(Mix, Nr, Nd);
Nr_t = Nr;
Nd_t = Nd;

%% Crop the signal
idxs = (crop_idxs+1) : 1 : (Nr - crop_idxs);
Txf = Txf_t(idxs,:);
Rxf = Rxf_t(idxs,:);
Tx  = Tx_t(idxs,:);
Rx  = Rx_t(idxs,:);
Mix = Mix_t(idxs,:);

[Nr, Nd] = size(Mix);



%% Plot frequencies
n_wave = 3;
figure(); hold on; box on; grid minor;
title(['Tx and Rx signal freq ',x0_str])
min_f = min(Txf_t(:)); max_f = max(Txf_t(:));
for i=1:n_wave
    xs = [((i-2)*Nr_t + idxs(end))*[1,1]*dt*1e9, ((i-1)*Nr_t + idxs(1))*[1,1]*dt*1e9];
    ys = [min_f,max_f,max_f,min_f];
    patch(xs,ys,'r','Facealpha',0.3,'EdgeColor', 'none');
    
    xs = [((i-1)*Nr_t + idxs(1))*[1,1]*dt*1e9, ((i-1)*Nr_t + idxs(end))*[1,1]*dt*1e9];
    ys = [min_f,max_f,max_f,min_f];
    patch(xs,ys,'green','Facealpha',0.3,'EdgeColor', 'none');
%     plot(((i-1)*Nr_t + idxs(1))*[1,1]*dt*1e9,[min_f,max_f],'--k')
%     plot(((i-1)*Nr_t + idxs(end))*[1,1]*dt*1e9,[min_f,max_f],'--k')
end
plot((1:Nr_t*n_wave)*dt*1e9, Txf_t(1:Nr_t*n_wave),'r','linewidth',3,'DisplayName','Tx')
plot((1:Nr_t*n_wave)*dt*1e9, Rxf_t(1:Nr_t*n_wave),'b','linewidth',3,'DisplayName','Rx')
xlabel("time [ns]");
ylabel("Freq [Hz]");
axis tight
set(gca,'FontSize',20)

%% 1D FFT to compute the range
% Compute fft and normalize
Y = fft(Mix, Nr, 1);
Y = Y/Nr;

% Take the absolute value of FFT output
P2 = abs(Y);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
P1 = P2(1:Nr/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:); % since we cutt half multiply x2 the non zero freq amplitude

% Compute the relative frequencies
Fs = 1/(dt);
f = Fs*(0:(Nr/2))/Nr;
d = speed_of_light*f/(2*slope);

%plotting the range
figure(); hold on; box on; grid minor;
title(['Range from 1D FFT ', x0_str])
plot(d',P1,'linewidth',2) 
% axis ([min(d) max(d) 0 0.5]);
xlabel("dist [m]");
ylabel("RxPwr [W]");
axis tight
set(gca,'FontSize',20)



%% RANGE DOPPLER RESPONSE
% This will run a 2DFFT on the mixed signal (beat signal) output and 
% generate a range doppler map.You will implement CFAR on the generated RDM

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(:,1:Nd);
sig_fft2 = fftshift (sig_fft2); % Keep bot side --> simmetry
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;


% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.
Fs = 1/(dt);
f = Fs*(-Nr/2:1:(Nr/2-1))/Nr;
range_axis = speed_of_light*f/(2*slope);

% Compute doppler_speed
Fs = 1/(dt*Nr_t);
f = Fs*(-Nd/2:1:(Nd/2-1))/Nd;
doppler_axis = f*wavelength/2;
% doppler_axis = linspace(-100,100,Nd);
% range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);

%% 2D Plot of RDM
figure, hold on;
title(['Range Doppler 2D FFT', x0_str]);
imagesc(doppler_axis,range_axis, RDM);
plot(x0(2),x0(1),'xr','markersize',25,'linewidth',3,'DisplayName','[s0,v0]')
c = colorbar;
c.Label.String = 'RxPwr [dBW]';
% surf(doppler_axis,range_axis,RDM, 'EdgeColor', 'none' );
xlabel("Doppler Velocity [m/s]");
ylabel("Range [m]");
legend
axis tight
set(gca,'FontSize',20)
%% CFAR 
%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 50;
Td = 5;

%Select the number of Guard Cells in both dimensions 
Gr = 10;
Gd = 2;

% offset the threshold by SNR value in dB
offset=pow2db(20);

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

% Vector to hold threshold values 
threshold_cfar = zeros(Nr,Nd);

%Vector to hold final signal after thresholding
signal_cfar = zeros(Nr,Nd);

% 2. Slide window across the signal length
for i = (Gr+Tr+1) : (Nr-(Gr+Tr))     
    for j = (Gd+Td+1) : (Nd-(Gd+Td))
        
        % 2. - 5. Determine the noise threshold by measuring it within the training cells
        % Calculate noise SUM in the area around CUT
        noise_level = 0;
        T = 0;
        for k = i-(Tr+Gr) : i+(Tr+Gr)
            for m = j-(Td+Gd) : j+(Td+Gd)
                if (abs(i-k) > Gr || abs(j-m) > Gd)
                    noise_level = noise_level + db2pow(RDM(k,m));
                    T = T + 1;
                end
            end
        end        
        threshold = offset + pow2db(noise_level/T);

        % 6. Measuring the signal within the CUT
        signal = RDM(i, j);

        % 8. Filter the signal above the threshold
        if signal<threshold
            signal = 0;
        else
            signal = 1;
        end
        threshold_cfar(i, j) = threshold;
        signal_cfar(i, j) = signal;
    end
end

%% 2D Plot filtered data
figure(); hold on; box; axis tight;
title(['Range Doppler Filtered 2D FFT', x0_str]);
imagesc(doppler_axis,range_axis,signal_cfar);
plot(x0(2),x0(1),'xr','markersize',20,'linewidth',2,'DisplayName','[s0,v0]');
plot([0,0],range_axis([1,end]),'--k','linewidth',3);
plot(doppler_axis([1,end]),[0,0],'--k','linewidth',3);
legend()
c = colorbar;
c.Label.String = 'Filtered RxPwr';
xlabel("speed [m/s]");
ylabel("distance [m]");
set(gca,'FontSize',20)

%% 3D plot of signal and threshold
figure;
title(['Range Doppler Threshold 2D FFT', x0_str]);
hold on; grid on;
surf(doppler_axis, range_axis, RDM, 'EdgeColor', 'none' );
tmp = threshold_cfar((Gr+Tr+1) : (Nr-(Gr+Tr)), :);
tmp = tmp(:, (Gd+Td+1) : (Nd-(Gd+Td)));
surf(doppler_axis((Gd+Td+1) : (Nd-(Gd+Td))), range_axis((Gr+Tr+1) : (Nr-(Gr+Tr))),...
    tmp, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none' );
legend("value","threshold");
% surf(doppler_axis,range_axis,RDM, 'EdgeColor', 'none' );
xlabel("Doppler Velocity [m/s]");
ylabel("Range [m]");
zlabel('RxPwr [dBW]');
set(gca,'FontSize',20)


%% FUNCTION
function [out, freq] = sawtoothWave(t,f0, Tchirp, Bsweep)
    f = 1/Tchirp;
    ks = 1:1000; % higher ks better approximation
    sign = 1; %-(2*rem(ks,2) - 1); % create alternatig signals
    freq = Bsweep*(-1/pi*sum(sign.*sin(2*pi*ks*f*t)./ks)+0.5) + f0;    
    phase = (Bsweep*(1/(2*pi*pi*f)*sum(sign.*cos(2*pi*ks*f*t)./(ks.*ks))+0.5*t) +f0*t)*2*pi;
    out = cos(phase);
end

 
 