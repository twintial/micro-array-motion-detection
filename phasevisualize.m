clc;clear;
close all;
DEBUG = 1;
% 构建声源
record = 1;
audioFrameLength = 2048;
if record
    audioFileName = 'audio/exercise/leg_raise1.wav';
    audioInput = dsp.AudioFileReader( ...
        'OutputDataType','double', ...
        'Filename',audioFileName, ...
        'PlayCount',inf, ...
        'SamplesPerFrame',audioFrameLength);
    fs = audioInput.SampleRate;
else
    fs = 48000;
    audioInput = audioDeviceReader(...
        'Device', 'miniDSP ASIO Driver', ...
        'Driver', 'ASIO', ...
        'SampleRate', fs, ...
        'NumChannels', 7 ,...
        'OutputDataType','double',...
        'SamplesPerFrame', audioFrameLength);
end

endTime = 50;
if ~record
    endTime = 1e4;
end

fc = 20e3;

dur = 15;
t = 1/fs:1/fs:dur;
cw_signal = cos(2*pi*fc*t); %生成余弦波信号，用来发射
cw_sin = sin(2*pi*fc*t);%解调时用的正弦波信号

% figure(1);
% y_g = [0 0];
% t_s = [0];
% graph_IQ = plot(t_s,y_g);
% title('v-t');
% xlabel('t(s)');
% ylabel('I/Q')
% NameArray = {'LineStyle'};
% ValueArray = {'-','--'}';
% set(graph_IQ,NameArray,ValueArray)

sgement_t = audioFrameLength/fs;
fc = 20e3;
% figure(2)
% h = animatedline('MaximumNumPoints', 10);

tic
for idx = 1:(endTime*fs/audioFrameLength)
    cycleStart = toc;
    y = audioInput();
    Wc = [2*(fc-3.5e3)/fs,2*(fc+3.5e3)/fs];
    [b, a] = butter(4,Wc);
    y = filter(b,a,y);
    y = y(:,1);
    
    % 转复信号
%     y2 = hilbert(y);
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
    data_signal = y1+1j.*y2;
    
%     ph = angle(data_signal);
%     ph = unwrap(ph);
%     figure(3);
%     t = 1/fs:1/fs:audioFrameLength/fs;
%     t_s = t + sgement_t*(idx-1);
%     plot(t_s, ph);

% %     圈圈图
%     real_y2 = real(y2);
%     imag_y2 = imag(y2);
%     addpoints(h, real_y2, imag_y2);
%     for i = 1:size(real_y2)
%         addpoints(h, real_y2(i), imag_y2(i));
%         drawnow
%         pause(0.01);
%     end
        
    real_y2 = y1;
    imag_y2 = y2;
    
    % real/imag图
    y_IQ = [real_y2;imag_y2];
%     y_g = [y_g; y_IQ];
    t = 1/fs:1/fs:audioFrameLength/fs;
%     t_s = [t_s t + sgement_t*(idx-1)];
    t_s = t + sgement_t*(idx-1);
    figure(1);
%     set(graph_IQ, 'XData', t_s, 'YData', y_g);
    plot(t_s,y_IQ)
    
    if record
        pause(audioFrameLength/fs - toc + cycleStart) % 为什么是这个值
    end
end
%% 离线，不分窗口
clc;clear;
close all;
audioFileName = 'audio/origin_mic_cross.wav';
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
y_1 = y;

% 转复信号
% data_signal = hilbert(y);
%解调IQ分量
for i=1:1:length(y)
    I_signal(i) = cw_signal(i)*y(i); %乘以余弦分量
    Q_signal(i) = cw_sin(i)*y(i);%乘以正弦信号
    
    I_signal_1(i) = cw_signal(i)*y_1(i); %乘以余弦分量
    Q_signal_1(i) = cw_sin(i)*y_1(i);%乘以正弦信号
end


%去掉IQ分量中的高频波，
Wn = 2*100/fs;
[b, a] = butter(4,Wn);
y1 = filter(b,a,I_signal);
y2 = filter(b,a,Q_signal);
data_signal = y1+1j.*y2;

% data_signal = data_signal./abs(data_signal) * 100;

y1_1 = filter(b,a,I_signal_1);
y2_1 = filter(b,a,Q_signal_1);
data_signal_1 = y1_1+1j.*y2_1;

ph = angle(data_signal);
ph = unwrap(ph);

ph_1 = angle(data_signal_1);
ph_1 = unwrap(ph_1);

% ph = ph - ph_1;
figure(1);
t = 1/fs:1/fs:length(data_signal)/fs;
plot(t(10000:end), ph(10000:end));


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
% g = plot(real_y2(1), imag_y2(1));
% set(gca, 'Xlim', [-0.025,0.025]);
% set(gca, 'Ylim', [-0.025,0.025]);
% t = 1/48000;
% for i = 48000*2:500:size(real_y2,2)
%     set(g, 'XData', real_y2(48e3*2:i), 'YData', imag_y2(48e3*2:i));
%     drawnow
%     disp(t*i)
% end


% h = animatedline('MaximumNumPoints', 1000);
% for i = 1:size(real_y2,2)
%     addpoints(h, real_y2(i), imag_y2(i));
%     drawnow
% end