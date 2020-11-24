%%
clc;
clear;
deviceReader1 = audioDeviceReader('Driver', 'DirectSound');
deviceReader2 = audioDeviceReader('Driver', 'ASIO');
devices_DirectSound = getAudioDevices(deviceReader1);
devices_ASIO = getAudioDevices(deviceReader2);
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
% deviceReader = audioDeviceReader();
% setup(deviceReader)

% 要加samplerate
fileWriter = dsp.AudioFileWriter('hand_move_face1.wav','FileFormat','WAV','SampleRate', fs);
disp('Speak into microphone now.');

tic;
while toc < 10
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