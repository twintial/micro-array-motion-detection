clear;clc;
close all;
%% Parameter Interface
Frequence0          = 60;        %单位：Hz    
Frequence1          = 130;       %单位：Hz
Frequence2          = 1e3;       %单位：Hz
SampleFre           = 4e3;       %单位：Hz
SampleLen           = SampleFre; %采样点数
%% Main
%-------------------产生三路信号
t = 0:1/SampleLen:1/SampleFre*(SampleLen-1);
SignalData0 = sin(2*pi*Frequence0*t);
SignalData1 = sin(2*pi*Frequence1*t);
SignalData2 = sin(2*pi*Frequence2*t);
SignalData3 = SignalData0+SignalData1+SignalData2;
figure;hold on
plot(t(1:150),SignalData0(1:150),'b')
plot(t(1:150),SignalData1(1:150),'r')
plot(t(1:150),SignalData2(1:150),'k')
hold off
figure;plot(t(1:150),SignalData3(1:150))

FFT_Data = fft(SignalData3);
L = SampleLen;
P2 = abs(FFT_Data/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = SampleFre*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% 高通滤波
% y = highpass(SignalData3,150,SampleFre);
% 滤波2
fs = SampleFre;
Wc = [2*100/fs,2*998/fs];   %end frequence 10hz
[b, a] = butter(8,Wc);
y = filter(b,a,SignalData3);

FFT_Data = fft(y);
L = SampleLen;
P2 = abs(FFT_Data/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = SampleFre*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
