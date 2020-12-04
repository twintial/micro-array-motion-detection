% 效果不好，可能是窗口太小
clc;clear;
audioFileName = 'audio/exercise/leg_raise1.wav';
audioFrameLength = 2048;
audioInput = dsp.AudioFileReader( ...
    'OutputDataType','double', ...
    'Filename',audioFileName, ...
    'PlayCount',inf, ...
    'SamplesPerFrame',audioFrameLength);
fs = audioInput.SampleRate;
fc = 20e3;

startTime = 12;
endTime = 40;
for i = 1:floor(startTime*fs/audioFrameLength)
    y = audioInput();
end
for idx = ceil(startTime*fs/audioFrameLength):(endTime*fs/audioFrameLength)
    y = audioInput();
    y(:,6) = [];
    % 滤波
    Wc = [2*(fc-3.5e3)/fs,2*(fc+3.5e3)/fs];
    [b, a] = butter(4,Wc);
    y = filter(b,a,y);
    y2 = hilbert(y);
    ph = angle(y2);
    ph = unwrap(ph);
    diff_ph = diff(ph);
    diff_ph = diff_ph(10:end-10,:); % 2028*6
    % 存储数据 
    temp = diff_ph';
    temp = temp(:)';
    save dataset/leg_raise1.txt temp -ascii -append
end

%% 不用hilbert
clc;clear;
close all;
audioFileName = 'audio/exercise/cross_streth_origin.wav';
[y, fs] = audioread(audioFileName);
fc = 20e3;
dur = 60;
t = 1/fs:1/fs:dur;
cw_signal = cos(2*pi*fc*t); %生成余弦波信号，用来发射
cw_sin = sin(2*pi*fc*t);%解调时用的正弦波信号

Wc = [2*(fc-3.5e3)/fs,2*(fc+3.5e3)/fs];
[b, a] = butter(4,Wc);
y = filter(b,a,y);

% y_1 = y(:,2);
y = y(:,1); % 取中心

% 转复信号
% data_signal = hilbert(y);
%解调IQ分量
for i=1:1:length(y)
    I_signal(i) = cw_signal(i)*y(i); %乘以余弦分量
    Q_signal(i) = cw_sin(i)*y(i);%乘以正弦信号
    
%     I_signal_1(i) = cw_signal(i)*y_1(i); %乘以余弦分量
%     Q_signal_1(i) = cw_sin(i)*y_1(i);%乘以正弦信号
end


%去掉IQ分量中的高频波，
Wn = 2*100/fs;
[b, a] = butter(4,Wn);
y1 = filter(b,a,I_signal);
y2 = filter(b,a,Q_signal);
data_signal = y1+1j.*y2;

% data_signal = data_signal./abs(data_signal);

% y1_1 = filter(b,a,I_signal_1);
% y2_1 = filter(b,a,Q_signal_1);
% data_signal_1 = y1_1+1j.*y2_1;

ph = angle(data_signal);
ph = unwrap(ph);

% ph_1 = angle(data_signal_1);
% ph_1 = unwrap(ph_1);

% ph = ph - ph_1;
figure(1);
t = 1/fs:1/fs:length(data_signal)/fs;
plot(t, ph);


real_y2 = real(data_signal);
imag_y2 = imag(data_signal);

% real/imag图
y_IQ = [real_y2;imag_y2];
%     y_g = [y_g; y_IQ];
% t = 1/fs:1/fs:audioFrameLength/fs;
%     t_s = [t_s t + sgement_t*(idx-1)];
t_s = t;
figure(2);
%     set(graph_IQ, 'XData', t_s, 'YData', y_g);
plot(t_s,y_IQ)


% 圈圈图
% figure(3)
% h = animatedline('MaximumNumPoints', 1000);
% for i = 1:size(real_y2,2)
%     addpoints(h, real_y2(i), imag_y2(i));
%     drawnow
% end
%% 文件读写操作
a = [1 2 3 4 5 6];
b = [7 8 9 10 11 12];
% fid = fopen('datatest.txt','w');
% fprint(fid, a, 'double');
% fclose(fid);
save datatest.txt a -ascii -append
save datatest.txt b -ascii -append