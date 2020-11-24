clc;clear;
close all;
fs = 48000; %����Ƶ��48KHz
dur = 10; %��������ʱ��10s
fc = 20000;%����Ƶ��ʱ17KHz
t = 1/fs:1/fs:dur;

cw_signal = cos(2*pi*fc*t); %�������Ҳ��źţ���������
cw_sin = sin(2*pi*fc*t);%���ʱ�õ����Ҳ��ź�
info = audiodevinfo;
player = audioplayer(cw_signal,fs,16,3);
% play(player);
% sound(cw_signal, fs);
R = audiorecorder(fs, 16,1,1) ;%����¼�ƶ��� ����¼�Ƶ��ǵ�����

play(player);
record(R,10); %%¼��ʱ��10s
pause(10);

stop(R);

pause(0.1);
stop(player);

myRecording = getaudiodata(R)';
result = myRecording; %¼��
y = result;
t1 = 1/fs:1/fs:length(y)/fs;
plot(t1,y);

%�˲�����ȡ�������ź�
Wc = [2*(fc-350)/fs,2*(fc+350)/fs];   %end frequence 10hz
[b, a] = butter(4,Wc);
y = filter(b,a,y);

%���IQ����
I_signal = ones(1);
Q_signal = ones(1);
for i=1:1:length(y)
    I_signal(i) = cw_signal(i)*y(i); %�������ҷ���
    Q_signal(i) = cw_sin(i)*y(i);%���������ź�
end
%ȥ��IQ�����еĸ�Ƶ����
Wn = 2*100/fs;
[b, a] = butter(4,Wn);
y1 = filter(b,a,I_signal);
y2 = filter(b,a,Q_signal);

%����λ�ͷ���
data_signal = y1+1j.*y2;
ph = angle(data_signal);
ph = unwrap(ph);
figure;
t1 = 1/fs:1/fs:length(y1)/fs;
plot(t1,ph);
