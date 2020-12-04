% music
clc;clear;
close all;

DEBUG = 0;
% 构建声源
record = 1;
audioFrameLength = 1024;
if record
    audioFileName = 'hand_move_face1.wav';
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
% 构建uca
microphone = phased.OmnidirectionalMicrophoneElement();
uca = phased.UCA('NumElements',6,'Radius',0.043);
if DEBUG
    viewArray(uca,'Title','Uniform Circular Array (UCA)','ShowIndex','All','ShowNormals',true);
end
% music
c = 343; % 声速
fc = 20e3; % 20kHz
amp = 100; % amplification
ns = 1; % 信号数量
musicazelspectrum = phased.MUSICEstimator2D('SensorArray',uca,...
    'OperatingFrequency',fc,'PropagationSpeed',c,...
    'AzimuthScanAngles',-90:90,'ElevationScanAngles',-1:30,...,
    'DOAOutputPort',true,...
    'NumSignalsSource','Property','NumSignals',ns);

endTime = 100;
DoaDisplayer = MyDoaDisplay(); % 可视化
% 循环处理信号
tic;
for idx = 1:(endTime*fs/audioFrameLength) % 为什么是这个值
    cycleStart = toc;
    y = audioInput();
    y = y(:,2:end);
    max_y = max(y);
    if min(max_y) == 0
        max_y = max_y + 1e-5;
    end
    y = y ./ max_y;
    % 画个小图
    t = 1/fs:1/fs:1024/fs;
    figure(1);
    for i = 1:size(y,2)
        y_i = y(:,i)';
        %     max_y = max(y_i);
        %     r = 1/max_y; % 放大倍率
        %     y_i = y_i * r;
        subplot(4,2,i);
        plot(t,y_i);
        xlabel('t');
        ylabel('A');
    end
    
    %滤波，提取超声波信号
    Wc = [2*(fc-350)/fs,2*(fc+350)/fs];
    [b, a] = butter(4,Wc);
    y = filter(b,a,y);
    % y = highpass(y,20e3-350,fs);
    [~,ang] = musicazelspectrum(y);
    figure(2);
    plotSpectrum(musicazelspectrum);
    drawnow
    %DoaDisplayer(deg2rad(ang(1))); % 可视化
    if record
        pause(audioFrameLength/fs - toc + cycleStart) % 为什么是这个值
    end
end
release(audioInput)

%% simu
clc;clear;
close all;

% 构建uca
microphone = phased.OmnidirectionalMicrophoneElement();
uca = phased.UCA('NumElements',6,'Radius',0.043);

ura = phased.URA('Size',[10 5],'ElementSpacing',[0.3 0.5]);

ang1 = [63; 45];         % First signal
ang2 = [-20; 20];        % Second signal

c = 343; % 声速
fc = 20e3; % 20kHz
% c = physconst('LightSpeed');
% fc = 300e6; % 20kHz
lambda = c / fc;
ns = 1; % 信号数量

Nsamp = 2048;
nPower = 0.01;
rs = rng(2007);

signal = sensorsig(getElementPosition(uca)/lambda,Nsamp, ...
    [ang1],nPower);
rng(rs);
% signal = signal ./ abs(max(signal));

musicazelspectrum = phased.MUSICEstimator2D('SensorArray',uca,...
    'OperatingFrequency',fc,'PropagationSpeed',c,...
    'AzimuthScanAngles',-90:90,'ElevationScanAngles',-10:60,...,
    'DOAOutputPort',true,...
    'NumSignalsSource','Property','NumSignals',ns);
[~,ang] = musicazelspectrum(signal);
plotSpectrum(musicazelspectrum);
disp(ang)

%% 波形图查看
clear;clc;
close all;
audioFileName = 'audio/test_face1.wav';
[y, fs] = audioread(audioFileName);
endTime = size(y,1)/fs;
t = 1/fs:1/fs:endTime;
y = y./max(y);
for i = 1:size(y,2)
    y_i = y(:,i)';
%     max_y = max(y_i);
%     r = 1/max_y; % 放大倍率
%     y_i = y_i * r;
    subplot(4,2,i);
    plot(t,y_i);
    xlabel('t');
    ylabel('A');
    y(:,i) = y_i';
end
% 保存
% fileWriter = dsp.AudioFileWriter('test_face1_after_amplify.wav','FileFormat','WAV','SampleRate', fs);
% step(fileWriter, y);
% release(fileWriter);

%% fft 能看到很好的现象。存在一个问题，如何判断周围是否有移动的物体。用相位，差分，存在问题，改用波峰。
clc;clear;
close all;

DEBUG = 1;
% 构建声源
record = 1;
audioFrameLength = 2048;
if record
    audioFileName = 'audio/exercise/curls3.wav';
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
c = 343; % 声速

% 定义图画对象
graph_fdB = graph;

v_s = [0];
t_s = [0];

% db阈值 16/10还不错，由于input6的问题这个值肯定也存在问题
threshold = 10;
% 10
thre_peak = 10;

if record == 0
    % 产生信号
    fs = 48000; %采样频率48KHz
    dur = 65; %发送声音时长10s
    fc = 20000;%中心频率时17KHz
    t = 1/fs:1/fs:dur;
    cw_signal = cos(2*pi*fc*t); %生成余弦波信号，用来发射
    info = audiodevinfo;
    player = audioplayer(cw_signal,fs,24,3);
    play(player);
end

disp("开始循环处理信号");
% 循环处理信号
tic;
for idx = 1:(endTime*fs/audioFrameLength) % 为什么是这个值
    cycleStart = toc;
    y = audioInput();
%     % power增强
%     max_y = max(y);
%     if min(max_y) == 0
%         max_y = max_y + 1e-5;
%     end
%     y = y ./ max_y;

    % 滤波
    Wc = [2*(fc-3.5e3)/fs,2*(fc+3.5e3)/fs];
    [b, a] = butter(4,Wc);
%     y = filter(b,a,y);
    y = y(:,1); % 只拿中心点
    % beamform
%     weight = ones(7, 1)/6;
%     weight(6) = 0;
%     y = ump8beamform(y, fs, weight);
    
    % fft
    FFT_Data = fft(y);
    L = audioFrameLength;
    P2 = abs(FFT_Data/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    dB = db(P1,'power'); % 转换为分贝
    
%     % 效果不是很好，明天试试速度的连续性
%     % 计算相位，用于看是否有运动的物体
%     y2 = hilbert(y);
%     ph = angle(y2);
%     ph = unwrap(ph);
%     diff_ph = diff(ph);
%     % 将前40个和后40个点拿去,平均住不管用
%     diff_ph = diff_ph(40:end-40);
%     diff_ph = smooth(diff_ph, 0.9); % 平滑化
%     diff_ph = diff_ph(50:end-50);
%     
%     diff_ph_std = std(diff_ph)*1e5; % 大于5就有速度，暂时先这样，不是很好
%     %disp(diff_ph_std);
%     if diff_ph_std < 6
%         addpoints(graph_vt,toc,0);
%         drawnow
%         continue;
%     end
    
    % 计算阈值中的max_dB
    range1 = (f >= fc-3.5e3) & (f <= fc+3.5e3);
    f1 = f(range1);
    dB1 = dB(range1);
    range2 = (f1 < fc-1e3) | (f1 > fc+1e3);
    f1 = f1(range2);
    dB1 = dB1(range2);
    max_dB = max(dB1);
    %disp(max_dB+threshold);
    
    % 计算fr
    range = (f >= fc-1e3) & (f <= fc+1e3);
    f_filter = f(range);
    dB_filter = dB(range);
    % fc_index = find(dB_filter == max(dB_filter)); % 20kHz的索引
    f_candidates = f_filter(dB_filter > (max_dB+threshold));
%     disp(max_dB+threshold);
    [~, fr_index] = max(abs(f_candidates - fc)); % 接收到的从手反射过来的频率
    fr = f_candidates(fr_index);
    % disp(fr);
    
    % 判断运动，这里的阈值不使用之前计算得到的阈值，由于input6的问题这个值肯定也存在问题
    p = findpeaks(dB_filter,'MinPeakHeight',max_dB+thre_peak);
    n_p = size(p,1);
    if n_p > 1
        % 有运动
        % 计算速度
        ft = fc;
        v = c * (fr-ft)/(fr+ft);
    else
        % 无运动
        v = 0;
    end
    v_s(end+1) = v;
    t_s(end+1) = audioFrameLength/fs * idx;
    
    if DEBUG
        % 画速度图
        if idx == 1
            figure(1);
            graph_vt = plot(t_s,v_s);
            title('v-t');
            xlabel('t(s)');
            ylabel('v(m/s)')
        else
            num = 100;
            if(size(v_s,2)>num)
                set(graph_vt, 'XData', t_s(end-num:end), 'YData', v_s(end-num:end));
            else
                set(graph_vt, 'XData', t_s, 'YData', v_s);
            end
            drawnow
        end
        
        % 画频谱图
        x = f_filter;
        y = dB_filter;
%         x = 1/fs:1/fs:(length(ph)-1)/fs;
%         y = diff_ph;
        if ~ishandle(graph_fdB)
            % graph = plot(f,P1);
            figure(2);
            graph_fdB = plot(x, y);
            title('Amplitude Spectrum of fr with 20kHz ft');
            xlabel('f (Hz)');
            ylabel('dB')
            
%             title('phase diff');
%             xlabel('t(s)');
%             ylabel('phase diff(rad)')
            % set(gca,'Xlim',[fc-3.5e3, fc+3.5e3]);
        else
            set(graph_fdB, 'XData', x, 'YData', y);
%             drawnow limitrate;
            drawnow
        end
%         figure(3)
%         findpeaks(dB_filter,x,'MinPeakHeight',max_dB+10);
    end
    if record
        pause(audioFrameLength/fs - toc + cycleStart) % 为什么是这个值
    end
end
disp("开始处理速度");
% 处理速度,只取10s-40s的
v_s_filter = v_s(t_s>10 & t_s<40);
% 计算一个速度用2048/48000s，划分成5s大小的窗口，0.2s为步长
for i = 10:0.2:35
    t_selected = t_s>=i & t_s<=i+5;
    v_window = v_s(t_selected);
    % 自相关
    [acf,lags] = autocorr(v_window, 'NumLags', 100);
end
release(audioInput)
%% beamform test
audioFrameLength = 2048;
audioFileName = 'test_face1.wav';
audioInput = dsp.AudioFileReader( ...
    'OutputDataType','double', ...
    'Filename',audioFileName, ...
    'PlayCount',inf, ...
    'SamplesPerFrame',audioFrameLength);
fs = audioInput.SampleRate;
afw = dsp.AudioFileWriter( ...
    'Filename', 'test_face1_beamformed.wav', ...
    'SampleRate', fs, ...
    'DataType', 'double');
endTime = 10;
for idx = 1:(endTime*fs/audioFrameLength)
    y = audioInput();
    weight = ones(7, 1);
    y = ump8beamform(y, fs, weight);
    afw(y);
end
release(audioInput);
release(afw);


% 输入一个实信号输出一个复信号，正交下变频
function [I, Q] = real2complex(sr, st, st_shift90)
I_origin = sr.*st;
Q_origin = sr.*st_shift90;
%低通滤波
Wn = 2*100/fs;
[b, a] = butter(4,Wn);
I = filter(b,a,I_origin);
Q = filter(b,a,Q_origin);
end

% ump-8的beamform,y.shape = (m,7),以中心点为基点合成
function y_b = ump8beamform(y, fs, weight)
% 固定角度，高度角为0，手掌正对1号mic
r = 0.043; % 43mm
theta = pi/6;
dist = [0, r, r*sin(theta), -r*sin(theta), -r, -r*sin(theta), r*sin(theta)];
delay = dist/343;
shifted_data = delayseq(y,delay,fs);
y_b = shifted_data * weight;
end