clc;clear;
close all;

DEBUG = 0;
% 构建声源
record = 0;
audioFrameLength = 1024;
if record
    audioFileName = 't.wav';
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
    'AzimuthScanAngles',-180:180,'ElevationScanAngles',-1:30,...,
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
    %滤波，提取超声波信号
    Wc = [2*(fc-350)/fs,2*(fc+350)/fs];
    [b, a] = butter(10,Wc);
    y = filter(b,a,y);
    % y = highpass(y,20e3-350,fs);
    [~,ang] = musicazelspectrum(y);
    plotSpectrum(musicazelspectrum);
    %DoaDisplayer(deg2rad(ang(1))); % 可视化
    if record
        pause(audioFrameLength/fs - toc + cycleStart) % 为什么是这个值
    end
end
release(audioInput)

%% simu

% 构建uca
microphone = phased.OmnidirectionalMicrophoneElement();
uca = phased.UCA('NumElements',6,'Radius',0.043);

ura = phased.URA('Size',[10 5],'ElementSpacing',[0.3 0.5]);

ang1 = [40; 45];         % First signal
ang2 = [-20; 20];        % Second signal

c = 343; % 声速
fc = 20e3; % 20kHz
% c = physconst('LightSpeed');
% fc = 300e6; % 20kHz
lambda = c / fc;
ns = 2; % 信号数量

Nsamp = 1024;
nPower = 0.01;
rs = rng(2007);

signal = sensorsig(getElementPosition(uca)/lambda,Nsamp, ...
    [ang1 ang2],nPower);
rng(rs);

musicazelspectrum = phased.MUSICEstimator2D('SensorArray',uca,...
    'OperatingFrequency',fc,'PropagationSpeed',c,...
    'AzimuthScanAngles',-90:90,'ElevationScanAngles',-10:60,...,
    'DOAOutputPort',true,...
    'NumSignalsSource','Property','NumSignals',ns);
[~,ang] = musicazelspectrum(signal);
plotSpectrum(musicazelspectrum);
disp(ang)

%% fft 能看到很好的现象。存在一个问题，如何判断周围是否有移动的物体。用相位，差分。
clc;clear;
close all;

DEBUG = 0;
% 构建声源
record = 1;
audioFrameLength = 2048;
if record
    audioFileName = 'hand_move.wav';
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
        'NumChannels', 8 ,...
        'OutputDataType','double',...
        'SamplesPerFrame', audioFrameLength);
end

endTime = 10;
if ~record
    endTime = 1e4;
end
fc = 20e3;
c = 343; % 声速

% 定义图画对象
graph_fdB = graph;

figure;
graph_vt = animatedline;
title('v-t');
xlabel('t');
ylabel('v')

% db阈值 16还不错
threshold = 15;
% 循环处理信号 idx=50开始出现波动
% for i = 1:50
%     y = audioInput();
% end
tic;
for idx = 1:(endTime*fs/audioFrameLength) % 为什么是这个值
    cycleStart = toc;
    y = audioInput();
    % 滤波
    Wc = [2*(fc-3.5e3)/fs,2*(fc+3.5e3)/fs];
    [b, a] = butter(4,Wc);
    y = filter(b,a,y);
    y = y(:,1); % 只拿中心点
    
    % fft
    FFT_Data = fft(y);
    L = audioFrameLength;
    P2 = abs(FFT_Data/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    dB = db(P1,'power'); % 转换为分贝
    
    % 效果不是很好，明天试试速度的连续性
    % 计算相位，用于看是否有运动的物体
    y2 = hilbert(y);
    ph = angle(y2);
    ph = unwrap(ph);
    diff_ph = diff(ph);
    % 将前40个和后40个点拿去,平均住不管用
    diff_ph = diff_ph(100:end-100);
    diff_ph = smooth(diff_ph, 0.5); % 平滑化
    diff_ph = diff_ph(50:end-50);
    
    diff_ph_std = std(diff_ph)*1e5; % 大于5就有速度，暂时先这样，不是很好
    %disp(diff_ph_std);
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
    
    % disp(max_dB+threshold);
    
    % 计算fr
    range = (f >= fc-1e3) & (f <= fc+1e3);
    f_filter = f(range);
    dB_filter = dB(range);
    
    % fc_index = find(dB_filter == max(dB_filter)); % 20kHz的索引
    f_candidates = f_filter(dB_filter > (max_dB+threshold));
    [~, fr_index] = max(abs(f_candidates - fc)); % 接收到的从手反射过来的频率
    fr = f_candidates(fr_index);
    
    % disp(fr);
    
    % 计算速度
    ft = fc;
    v = c * (fr-ft)/(fr+ft);
    
    % 画v-t图
    addpoints(graph_vt,toc,v)
    % drawnow limitrate;
    drawnow
    
    if DEBUG
        x = 1/fs:1/fs:(length(ph)-299)/fs;
        %         y = dB_filter;
        y = diff_ph;
        if ~ishandle(graph_fdB)
            % graph = plot(f,P1);
            figure;
            graph_fdB = plot(x,y);
            title('Single-Sided Amplitude Spectrum of X(t)');
            xlabel('f (Hz)');
            % ylabel('|P1(f)|');
            ylabel('dB')
            % set(gca,'Xlim',[fc-3.5e3, fc+3.5e3]);
        else
            set(graph_fdB, 'XData', x, 'YData', y);
            % drawnow limitrate;
            drawnow
        end
    end
    if record
        pause(audioFrameLength/fs - toc + cycleStart) % 为什么是这个值
    end
    pause(0.001);
end

% 输入一个实信号输出一个复信号
function [I, Q] = real2complex(sr, st, st_shift90)
I_origin = sr.*st;
Q_origin = sr.*st_shift90;
%低通滤波
Wn = 2*100/fs;
[b, a] = butter(4,Wn);
I = filter(b,a,I_origin);
Q = filter(b,a,Q_origin);
end