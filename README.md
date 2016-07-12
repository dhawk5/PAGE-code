%% TEST intensity_estimates_func() (ief)
% This code generates simulated data from a plane wave propagating in the
% +y direction and feeds the data into the intensity_estimates_func() to
% see if the correct values for the equilateral triangle probe are achieved
% for the various energy quantities.
% Modified July 2016 for 1d 3 mic probe - TBN
clear; 
close all;
addpath 'P:\MATLAB\Code Directory'
%myFigureDefaults(2,'lg') % For papers Set text scaling to 1 and text to 'lg' style

% Define Constants
c = 343;
rho = 1.21;
Iref = 1e-12;
Pref = 20e-6;
%% Select channels, ID num and spacing
dataDir = 'P:\MATLAB\Insulated Driver 6-30 and 7-6';
ch = [0,3];%2,3]; % channel numbers
IDnum = 2; % ID number
IDstr = 'ID';
angle = 120;
spacing = 12 % in inches

% Position of probe
ppos = [0,1,0]; % (m)
% distance from
a = 0.0254*spacing; % in m
% Probe geometry
if length(ch) == 2
probe_config = [ 0 0 0
    0 a 0];
elseif length(ch) == 3
probe_config = [0 -a 0 
    0 0 0
    0 a 0];
end
% Equilateral trangle 2d probe
% probe_config = a*[0,0,0;...
%     0,-1,0;...
%     -sqrt(3)/2,1/2,0;...
%     sqrt(3)/2,1/2,0];




rvec = repmat(ppos,size(probe_config,1),1)+probe_config;

% We only need the y component because it is a plane wave
y_mic = rvec(:,2);

% Set up the wavenumber k
fs = 50000;
ns = 25000;
df = fs/ns;
% fss = 0:df:(fs/2-df);


w = hann(ns);
W = mean(w.*conj(w));   % Used to scale the ASD for energy conservation



% One-dimensional probe 
% deltaR = r(ch(3)+1,:) - r(ch(1)+1,:);
% R = deltaR(2); % For a 1-d probe along the y axis
for jch = 1:length(ch)
    x = binfileload(dataDir,IDstr,IDnum,ch(jch));
    N = length(x);
    [Xss(jch,:,:),numblocks,fss] = computeBlockFFt_dh(x,ns,w,W,fs);
    % Xss has dimensions numchannels x numfreqs x numblocks  - TBN
    % When time-averaging, the sum is over blocks NOT freqs
end
k = 2*pi*fss/c;

Xss = permute(Xss,[1 3 2]);

% Check the spectral levels - TBN
XssdB = 20*log10(squeeze(rms(Xss,2)/Pref));
figure()
semilogx(fss,XssdB(1,:))
hold on
semilogx(fss,XssdB(2,:))
if length(ch) >=3
semilogx(fss,XssdB(3,:))
end
%semilogx(f,20*log10(rms(Xss(4,:,:),3)/pref))
legend('1','2','3')
xlabel('Frequency (Hz)')
ylabel('SPL (dBre 20 \muPa')
title({'Insulated Driver'})
xlim([min(fss) 2e4])
ylim([-35 105]) 
grid on

%% Test intensity_estimates_func
[TRAD,PAGE] = intensity_estimates_func(fss,Xss,probe_config,rho,c,1);

%Analytical Plane Wave
P2 = squeeze(abs(Xss(2,:,:)));
Prms = rms(P2,1); % mean of the average pressures over blocks for center microphone, one value for each frequency - TBN
IaAn = Prms.^2/rho/c;
IaAndB = 10*log10(IaAn/Iref);
%% % plot intensities - I is the active intensity and Q is the reactive intensity
% Is abs needed? Appear to be some negative I values...
figure
PAGE.IdB = 10*log10((PAGE.Imag)/Iref);
TRAD.IdB = 10*log10((TRAD.Imag)/Iref);
semilogx(fss,IaAndB)
hold on
semilogx(fss,TRAD.IdB(1,:))
semilogx(fss,PAGE.IdB(1,:))
xlabel('Frequency (Hz)')
ylabel('Intensity (dB re 1pW/m^2)')
xlim([min(fss) 4e3])
ylim([0 100]) 
grid on
legend('Analytical','Traditional','PAGE')
title('Active Intensity Calculations')
filename = ['QM2Intensity',num2str(angle),'deg',num2str(length(ch)),'mic'];
saveFigure([dataDir,'/',filename])
%% plot direction
figure()
semilogx(fss,PAGE.Idir(1,:))
hold on
semilogx(fss,TRAD.Idir(1,:))
xlabel('freq (Hz)')
ylabel('Phase Difference (\circ)')
xlim([20 4e3])
ylim([-10 10])
legend('PAGE','Traditional')
ax = gca;
ax.YTick = [-10 -5 0 5 10];
grid on
title('Direction')
%% 



% %% Using a plane wave, create data to send into the ief
% % Plane wave propagating in the +y direction $p=A*exp(-1j*k*y)$
% 
% 
% % Amplitude is 1 Pa, so SPL should be ~91 dB
% A = 1;
% 
% % Calculate pressures in frequency domain
% Xss = permute(A*exp(-1j*k'*y_mic')',[1,3,2])/sqrt(2);
% % Since Gxxdf had to be divided by 2, we better divide Xss by sqrt(2).
% % Since we want conj(Xss).*Xss to equal Gxxdf in the
% % intensity_estimates_func()
% 
% %% Determine analytical values for Energy quantities
% 
% % SPL
% SPL_exact = 20*log10(A/(sqrt(2)*Pref));
% 
% % analytic particle velocity
% v = A/(rho*c)*exp(-1j*k'*y_mic')';
% 
% % Analytic active intensity
% Ia = A^2/(2*rho*c);
% SILa_exact = 10*log10(Ia/Iref)*ones(size(fss));
% 
% % Analytic reactive intensity
% Ir = zeros(size(fss));
% 
% % Analytic Potential Energy
% Ep = abs(A)^2/(4*rho*c^2);
% 
% % Analytic Kinetic Energy
% Ek = Ep;
% 
% % Analytic Impedance
% z = rho*c;
% 
% %% Expected Error
% % SPL no error for either method-Therefore Ep no error for either method
% % Ir should be zero for both methods
% % Ia, z, and Ek should be exact for the PAGE and fall off for the FD
% 
% %% Test SPL
% Gxxdf = conj(Xss).*Xss; % Pa^2
% SPL = 10*log10(squeeze(Gxxdf(1,:,:))/4e-10);
% assert(all(abs(SPL-SPL_exact)<0.2),'SPL is wrong')
% 
% %% Test intensity_estimates_func
% [FD,PAGE] = intensity_estimates_func(fss,Xss,probe_config,rho,c,1);
% 
% % Test PAGE intensity
% figure()
% plot(fss, PAGE.Imag,fss,Ia*ones(size(fss)),'--')
% xlabel('Freq (Hz)')
% ylabel('Intensity')
% legend('PAGE','Analytical')
% if all(10*log10(PAGE.Imag./Ia)<0.2)
%     disp('SIL Active PAGE passed')
% else
%     disp('SIL Active PAGE is wrong')
% end
% 
% % Test PAGE reactive intensity
% if any(any(PAGE.Q))
%     disp('SIL Reactive PAGE passed')
% else
%     disp('SIL Reactive PAGE is wrong')
% end
% 
% % Test PAGE Potential Energy
% if all(10*log10(PAGE.Ep/Ep)<0.2)
%     disp('Ep PAGE passed')
% else
%     disp('Ep PAGE is wrong')
% end
% 
% % Test PAGE Kinetic Energy
% if all(10*log10(PAGE.Ek/Ek)<0.2)
%     disp('Ek PAGE passed')
% else
%     disp('Ek PAGE is wrong')
% end
% 
% % Test Traditional intensity
% % figure
% % plot(fss, FD.Imag*2, fss,Ia*ones(size(fss)),fss,abs(Ia*2.*(sin(k*a/2)+sin(k*a))./(3*k*a)),'--')
% % xlabel('Freq (Hz)')
% % ylabel('Intensity')
% % legend('ief','Analytical','Expected error')
% if all(10*log10(FD.Imag(2:end)./(Ia*2*(sin(k(2:end)*a/2)+sin(k(2:end)*a))./(3*k(2:end).*a)))<0.2)
%     disp('SIL Active Traditional passed')
% else
%     disp('SIL Active Traditional is wrong')
% end
% 
% % Test Traditional reactive intensity
% if any(any(FD.Q))
%     disp('SIL Reactive Traditional passed')
% else
%     disp('SIL Reactive Traditional is wrong')
% end
% 
% % Test Traditional Potential Energy
% if all(10*log10(FD.Ep/Ep)<0.2)
%     disp('Ep Traditional passed')
% else
%     disp('Ep Traditional is wrong')
% end
% 
% % Test Traditional Kinetic Energy
% figure()
% plot(fss,FD.Ek,fss,Ek*ones(size(fss)),fss,Ek*sinc(3*k*a/(4*pi)).^2,'--')
% xlabel('Frequency (Hz)')
% ylabel('Kinetic Energy')
% legend('Traditional','Analytical','Expected Traditional')
% if all(10*log10(FD.Ek(2:end)/(Ek*sinc(3*k(2:end)*a/(4*pi)).^2))<0.2)
%     disp('Ek Traditional passed')
% else
%     disp('Ek Traditional is wrong')
% end
% 
% 
% 
% 
% 
% 
% 
% 
