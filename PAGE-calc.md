%% Load in data from plane-wave tube measurements
clc
close all;
clear

%define constants
path = 'P:\MATLAB\Plane Wave Tube 6-09';
fs = 48000;
ns = 2^13;
df = fs/ns;
fs1 = 0:df:(fs/2-df);
w = hann(ns);
W = mean(w.*conj(w));   % Used to scale the ASD for energy conservation
ch = [1,2,3,4];         % Which two mics will be used
r1 = [0 0 0];           % Position of first mic
r2 = [0 2 0];           % Position of second mic
r3 = [0 6 0];           % Position of third mic
r4 = [0 12 0];          % Position of fourth mic
rho = 1.21;             % Air density
c = 343;                % Speed of sound
Iref = 1e-12;           % Reference intensity
pref = 2e-5;            % Reference pressure
omega = 2*pi;           % Frequency array%

%x1 = binfileload(path,'ID',3,0);
%x2 = binfileload(path,'ID',3,1);
%x3 = binfileload(path,'ID',3,2);
%x4 = binfileload(path,'ID',3,3);

%  loglog(f,coh)
%  ylim([.1 1])
% When using binfileload:
% 3 - Anechoic broadband
% 4 - Anechoic filtered
% 5 - Rigid broadband
% 6 - Rigid filtered

deltaR = r2 - r1;
R = deltaR';

for ch = [1,2]
    x = binfileload(path,'ID',5,ch-1);
    N = length(x);
    [Xss(ch,:,:),numblocks,f] = computeBlockFFt_dh(x,ns,w,W,fs);
    % Xss has dimensions numchannels x numfreqs x numblocks  - TBN
    % When time-averaging, the sum is over blocks NOT freqs
end

% Check the spectral levels - TBN
figure()
semilogx(f,20*log10(rms(Xss(1,:,:),3)/pref))
hold on
semilogx(f,20*log10(rms(Xss(2,:,:),3)/pref))
legend('1', '2')
xlabel('Frequency (Hz)')
ylabel('SPL (dBre 20 \muPa')
xlim([20 3e3])

%initializing r-vector
% mic_config = [0,0,0;0,2,0;0,6,0;0,12,0];
% rvec = repmat(ppos,size(mic_config,1),1)+mic_config;
% y_mic = rvec(:,2);
% X1 = [y_mic(2)-y_mic(1);y_mic(3)-y_mic(2);y_mic(4)-y_mic(3);...
%    y_mic(4)-y_mic(2);y_mic(4)-y_mic(1)];

%% Calculate Avg. Pressure and Phase
P1 = squeeze(abs(Xss(1,:,:)));%rms(Xss(1,:,:),3); 
%P1 = squeeze(P1);
P2 = squeeze(abs(Xss(2,:,:)));
%P2 = squeeze(P2);
Pavg = mean((P1 + P2)/2,2); % mean of the average pressures over blocks, one value for each frequency - TBN
deltaP = (mean((P2 - P1),2));% mean of the gradient over blocks, one value for each frequency - TBN


Phi1 = atan2(imag(Xss(1,:,:)),real(Xss(1,:,:)));
Phi2 = atan2(imag(Xss(2,:,:)),real(Xss(2,:,:)));
deltaPhi = squeeze(Phi2 - Phi1);
deltaPhiAvg = mean(deltaPhi,2); % mean of deltaPhi over blocks, one value for each frequency - TBN


%% Use calculations to find Intensities
%active intensity Ia
 gradPhi = deltaPhiAvg/R(2);% simpler form for two mic probe (1/(R(2)*R(2)))*R(2) * PhiAvg; 
 Psquared = (Pavg.^2);
 omega =  2*pi*f';
 Ia = 1/rho./omega.*Psquared.*gradPhi;  % change f to omega - TBN
 Ia = 10*log10(Ia/Iref);
 
 
%reactive intensity Ir
gradP = deltaP/R(2); % simpler for two mic probe (1/(R(2)*R(2)))*R(2) * deltaP;
Ir = -1/rho./omega.*Pavg.*gradP; % change f to omega - TBN
Ir = 10*log10(Ir/Iref);

% plot intensities
figure()
semilogx(f,abs(Ia))
%title('Active Intensity (I_a)')
xlabel('freq (Hz)')
ylabel('Intensity (dB re 1pW/m^2)')
xlim([20 4e3])
ylim([20 70])

hold on
%figure()
semilogx(f,abs(Ir))
%title('Reactive Intensity (I_r)')
legend('I_a','I_r')
xlabel('freq (Hz)')
ylabel('Intensity (dB re 1pW/m^2)')
xlim([20 4e3])
ylim([20 70])
%% Plot Measurements of different microphones %
% [Xss1,numblocks,f] = computeBlockFFt_dh(x1,ns,w,W,fs);
% [Xss2,numblocks,f] = computeBlockFFt_dh(x2,ns,w,W,fs);
% [Xss3,numblocks,f] = computeBlockFFt_dh(x3,ns,w,W,fs);
% [Xss4,numblocks,f] = computeBlockFFt_dh(x4,ns,w,W,fs);
%  [f,G11] = autospec(x1,fs, ns, N, 1);
%  [f,G22] = autospec(x2,fs, ns, N, 1);
% [f,G33] = autospec(x3,fs, ns, N, 1);
% [f,G44] = autospec(x4,fs, ns, N, 1);

%% Checking coherence
% x1 = binfileload(path,'ID',3,0);
% x2 = binfileload(path,'ID',3,1);
% [f,G11] = autospec(x1,fs, ns, N, 1);
%  [f,G22] = autospec(x2,fs, ns, N, 1);
%  [f,G12] = crossspec(x1,x2,fs,ns,N,1);
%  coh = abs(G12).^2./(G11.*G22);
%  loglog(f,coh)
%  ylim([.1 1])

% figure
% semilogx(f,10 * log10(G11/4e-10),'r')
% hold on
% grid on
% title('first microphone')
% xlim([20 2e4])
% xlabel('freq (Hz)')
% ylabel('Autospectra (dB-re 20\muPa)')
% semilogx(f,10 * log10(G22/4e-10),'b')
% semilogx(f,10 * log10(G33/4e-10),'g')
% semilogx(f,10 * log10(G44/4e-10),'k')
% legend('0 in','2 in','6 in','12 in')
