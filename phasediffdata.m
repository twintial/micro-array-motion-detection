clc;clear;
audioFileName = 'audio/exercise/shoulder_rotation1.wav';
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
    save dataset/shoulder_rotation1.txt temp -ascii -append
end

%% 文件读写操作
a = [1 2 3 4 5 6];
b = [7 8 9 10 11 12];
% fid = fopen('datatest.txt','w');
% fprint(fid, a, 'double');
% fclose(fid);
save datatest.txt a -ascii -append
save datatest.txt b -ascii -append