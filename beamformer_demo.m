clear
clc;

DEBUG = 0;
REC = 0;
TEST_LENGTH = 30; % Streaming duration, in sec


%% ### Global parameters and devices ###
fs = 16000;
amp = 100; % amplification
c = 343; % sound speed, in m/s
audioFrameLength = 1024;
deviceReader = audioDeviceReader(...
 'Device', 'miniDSP ASIO Driver',...
 'Driver', 'ASIO', ...
 'SampleRate', fs, ...
 'NumChannels', 16 ,...
 'OutputDataType','double',...
 'SamplesPerFrame', audioFrameLength);

% Remap the sources to match Matlab's URA indexing
deviceReader.ChannelMappingSource = 'Property';
deviceReader.ChannelMapping = [8 6 4 2 7 5 3 1 10 12 14 16 9 11 13 15];
setup(deviceReader);

if REC
    % set up audio device writer
    audioWriter = audioDeviceWriter('SampleRate',fs, ...
            'SupportVariableSizeInput', false);
    isAudioSupported = (length(getAudioDevices(audioWriter))>1);
    fileWriter = dsp.AudioFileWriter(...
        'UMA16_rec.wav',...
        'FileFormat','WAV');
end


%% ### UMA-16 definition ##
M = 4;   % Number of elements on each row
N = 4;   % Number of elements on each column
dy = 0.042; % Spacing between elements on each row (m)
dz = 0.042; % Spacing between elements on each column (m)

microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',...
    [20 20e3]);
ura = phased.URA('Element',microphone,'Size',[N M],...
    'ElementSpacing',[dz dy]);

if DEBUG
    viewArray(ura,'Title','Uniform Rectangular Array (URA)');
end

% ### Beamformer setup ###
angSteer = [30;0]; % Steering angles in deg. (azimuth, elevation)
beamformer = ...
    phased.GSCBeamformer('SensorArray',ura,'SampleRate',fs,...
    'PropagationSpeed',c,'FilterLength',20,'Direction', angSteer);


%%
if DEBUG
    % Sinks for the streaming signals
    mic1_log = dsp.SignalSink;
    bf_log = dsp.SignalSink;
end

disp('Beginning streaming.')
tic;
while toc < TEST_LENGTH
    acq = amp * deviceReader();
    temp = beamformer(acq);
    
    play(audioWriter, temp);
    
    if REC
        fileWriter(acq);
    end
    
    if DEBUG
        bf_log(temp);
        mic1_log(acq(:,1));
    end
end
disp('Streaming finished.')

release(deviceReader);
if REC
    release(audioWriter);
end
