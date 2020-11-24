clc;clear;
close all;
fs = 48000; %采样频率48KHz
dur = 10; %发送声音时长10s
fc = 20000;%中心频率时17KHz
t = 1/fs:1/fs:dur;

cw_signal = cos(2*pi*fc*t); %生成余弦波信号，用来发射
cw_sin = sin(2*pi*fc*t);%解调时用的正弦波信号
info = audiodevinfo;
player = audioplayer(cw_signal,fs,16,3);
% play(player);
% sound(cw_signal, fs);
R = audiorecorder(fs, 16,1,1) ;%创建录制对象 这里录制的是单声道

play(player);
record(R,10); %%录音时长10s
pause(10);

stop(R);

pause(0.1);
stop(player);

myRecording = getaudiodata(R)';
result = myRecording; %录音
y = result;
t1 = 1/fs:1/fs:length(y)/fs;
plot(t1,y);

%滤波，提取超声波信号
Wc = [2*(fc-350)/fs,2*(fc+350)/fs];   %end frequence 10hz
[b, a] = butter(4,Wc);
y = filter(b,a,y);

%解调IQ分量
I_signal = ones(1);
Q_signal = ones(1);
for i=1:1:length(y)
    I_signal(i) = cw_signal(i)*y(i); %乘以余弦分量
    Q_signal(i) = cw_sin(i)*y(i);%乘以正弦信号
end
%去掉IQ分量中的高频波，
Wn = 2*100/fs;
[b, a] = butter(4,Wn);
y1 = filter(b,a,I_signal);
y2 = filter(b,a,Q_signal);

%求相位和幅度
data_signal = y1+1j.*y2;
ph = angle(data_signal);
ph = unwrap(ph);
figure;
t1 = 1/fs:1/fs:length(y1)/fs;
plot(t1,ph);
