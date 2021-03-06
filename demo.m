%%
clc;
clear;
deviceReader1 = audioDeviceReader('Driver', 'DirectSound');
deviceReader2 = audioDeviceReader('Driver', 'ASIO');
deviceWriter = audioDeviceWriter('Driver', 'ASIO');
devices_DirectSound = getAudioDevices(deviceReader1);
devices_ASIO = getAudioDevices(deviceReader2);
devices_output = getAudioDevices(deviceWriter);
info = audiodevinfo;
% 
% playRec = audioPlayerRecorder;
% asiosettings()

%% 录音
clc;
clear;
fs = 48000;
audioFrameLength = 2048;
deviceReader = audioDeviceReader(...
    'Device', 'miniDSP ASIO Driver', ...
    'Driver', 'ASIO', ...
    'SampleRate', fs, ...
    'NumChannels', 7 ,...
    'OutputDataType','double',...
    'SamplesPerFrame', audioFrameLength);
% deviceReader = audioDeviceReader(...
%     'SampleRate', fs, ...
%     'SamplesPerFrame', audioFrameLength);
% setup(deviceReader)

% 要加samplerate
fileWriter = dsp.AudioFileWriter('audio/daily_test.wav','FileFormat','WAV','SampleRate', fs);
disp('Speak into microphone now.');

% 产生信号
fs = 48000; %采样频率48KHz
dur = 10; %发送声音时长10s
fc = 20000;%中心频率时17KHz
t = 1/fs:1/fs:dur;
cw_signal = cos(2*pi*fc*t); %生成余弦波信号，用来发射
info = audiodevinfo;
player = audioplayer(cw_signal,fs,24,3);
play(player)

tic;
while toc < dur+5
    disp(toc)
    acquiredAudio = record(deviceReader);
    step(fileWriter, acquiredAudio);
    % fileWriter(acquiredAudio);
end
disp('Recording complete.');

release(deviceReader)
release(fileWriter)

%%
% beamforming
clc;clear;
close all;
audioFileName = 'mySpeech.wav';
audioFrameLength = 1024;
audioInput = dsp.AudioFileReader( ...
    'OutputDataType','double', ...
    'Filename',audioFileName, ...
    'PlayCount',inf, ...
    'SamplesPerFrame',audioFrameLength);
fs = audioInput.SampleRate;
fileWriter = dsp.AudioFileWriter('beamform.wav','FileFormat','WAV');
endTime = 5;
for idx = 1:(endTime*fs/audioFrameLength)
%     % 单纯累加
%     multichannelAudioFrame = audioInput();
%     x = sum(multichannelAudioFrame,2);
%     step(fileWriter, x);
end
release(audioInput);
release(fileWriter);
%%
clear;
clc;
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
    fs = 44100;
    audioInput = audioDeviceReader(...
        'Device', 'miniDSP ASIO Driver', ...
        'Driver', 'ASIO', ...
        'SampleRate', fs, ...
        'NumChannels', 8 ,...
        'OutputDataType','double',...
        'SamplesPerFrame', audioFrameLength);
end
% mic array
ump_8 = ump_8_sp_2d(6, 0.043, 'ump8');
% wavelength
c = 343; % sound speed, in m/s
f = 2e4;
wavelength = c / f;

endTime = 200;
DoaDisplayer = MyDoaDisplay();
tic;
for idx = 1:(endTime*fs/audioFrameLength) % 为什么是这个值
%     if idx == 5
%         break
%     end
    cycleStart = toc;
    multichannelAudioFrame = audioInput();
    multichannelAudioFrame = highpass(multichannelAudioFrame,2e4-2e3,fs);
    X = multichannelAudioFrame';
    X = X(1:7,:);
    R = (X*X')/ audioFrameLength;
    % normal MUSIC
    sp_normal = music_1d(R, 1, ump_8, wavelength, 1800, 'RefineEstimates', true);
    %sp_normal.true_positions = linspace(-pi/8, pi/3, 8);
    fprintf('[MUSIC] Estimated DOAs:\n');
    disp(sp_normal.x_est);
    DoaDisplayer(sp_normal.x_est);
    if record
        pause(audioFrameLength/fs - toc + cycleStart) % 为什么是这个值
    end
    %plot_sp(sp_normal, 'title', 'Normal MUSIC', 'PlotType', 'Polar');
end
release(audioInput)
%%

% design and visualize arrays
clc
clear(); close all;
% 
wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 2, d);
% design_cp = design_array_1d('coprime', [4 5], d);
% design_nested = design_array_1d('nested', [5 7], d);
% design_mra = design_array_1d('mra', 12, d);
% 
visualize_array(design_ula, 'VisualizeCoarray', true);
% visualize_array(design_cp, 'VisualizeCoarray', true);
% visualize_array(design_nested, 'VisualizeCoarray', true);
% visualize_array(design_mra, 'VisualizeCoarray', true);

design_uca = uca_2d(6, 10, 'uca');
visualize_array(design_uca, 'VisualizeCoarray', true);
% 
% ump_8 = ump_8_sp_2d(6, 10, 'ump8');
% visualize_array(ump_8, 'VisualizeCoarray', true);

%%
% using MDL and AIC to find source number, and MUSIC to find DOAs
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 12, d);
doas = linspace(-pi/8, pi/3, 8);
power_source = 1;
power_noise = 1;
snapshot_count = 1000;
source_count = length(doas);

% stochastic (unconditional) model
[~, R] = snapshot_gen_sto(design_ula, doas, wavelength, snapshot_count, power_noise, power_source);

% source number detection
[~, l] = eig(0.5*(R+R'), 'vector');
l = flipud(l);
n_mdl = sn_mdl(l, design_ula.element_count, snapshot_count);
n_aic = sn_aic(l, design_ula.element_count, snapshot_count);
fprintf('There are %d sources.\n', source_count);
fprintf('# of sources estimated by MDL = %d\n', n_mdl);
fprintf('# of sources estimated by AIC = %d\n', n_aic);

% normal MUSIC
tic;
sp_normal = music_1d(R, source_count, design_ula, wavelength, 180, 'RefineEstimates', true);
toc;
sp_normal.true_positions = doas;
fprintf('[MUSIC] Estimated DOAs:\n');
disp(sp_normal.x_est);
plot_sp(sp_normal, 'title', 'Normal MUSIC', 'PlotType', 'Polar');

% root MUSIC
tic;
sp_root = rmusic_1d(R, source_count, 2*pi*design_ula.element_spacing/wavelength);
toc;
sp_root.true_positions = doas;
fprintf('[Root-MUSIC] Estimated DOAs:\n');
disp(sp_root.x_est);
plot_sp(sp_root, 'Title', 'Root-MUSIC');


%%
clc; clear; close all;
%%%%%%%%%%%%%%%%%% MUSIC %%%%%%%%%%%%%%%%%%
derad = pi/180;
N = 8;               % 阵元个数        
M = 3;               % 信源数目
theta = [-50 0 50];  % 待估计角度
snr = 10;            % 信噪比
K = 1024;            % 快拍数 10/1024

dd = 0.5;            % 阵元间距 d=lamda/2      
d=0:dd:(N-1)*dd;
A=exp(-1i*2*pi*d.'*sin(theta*derad));

S=randn(M,K); X=A*S;
X1=awgn(X,snr,'measured');
Rxx=X1*X1'/K;
%InvS=inv(Rxx); 
[EV,D]=eig(Rxx);    % 特征向量 特征值
EVA=diag(D)';
[EVA,I]=sort(EVA);  % 从小到大排列 返回索引
%EVA=fliplr(EVA);   % 反转元素
EV=fliplr(EV(:,I)); % 按列反转特征向量

for iang = 1:361    % 遍历
        angle(iang)=(iang-181)/2;
        phim=derad*angle(iang);
        a=exp(-1i*2*pi*d*sin(phim)).';
        L=3;    
        En=EV(:,L+1:N);
        %SP(iang)=(a'*a)/(a'*En*En'*a);
        SP(iang)=1/(a'*En*En'*a);
end
SP=abs(SP);
%SPmax=max(SP);
%SP=10*log10(SP/SPmax);  % 归一化
SP=10*log10(SP);
h=plot(angle,SP);
set(h,'Linewidth',0.5);
xlabel('入射角/(degree)');
ylabel('空间谱/(dB)');
%axis([-100 100 -40 60]);
set(gca, 'XTick',[-100:20:100]);
grid on;  

%% 自相关测试
clc;clear;
close all;
fs = 10;
x = 0:1/fs:200/fs;
y = sin(x);
plot(y);
figure(2);
[acf,lags,bounds,h] = autocorr(y, 'NumLags', 199);
% figure(3);
% [c,lags] = xcorr(y,y, 199);
% stem(lags,c);