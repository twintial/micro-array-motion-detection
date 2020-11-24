clf;

fs=100;N=256;   %采样频率和数据点数
n=0:N-1;t=n/fs;   %时间序列
x=0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t); %信号

y1=fft(x,N);    %对信号进行快速Fourier变换
y2=fftshift(y1);

mag1=abs(y1);     %求得Fourier变换后的振幅
mag2=abs(y2);    

f1=n*fs/N;    %频率序列
f2=n*fs/N-fs/2;

subplot(3,1,1),plot(f1,mag1,'r');  %绘出随频率变化的振幅
xlabel('频率/Hz');
ylabel('振幅');title('图1：usual FFT','color','r');grid on;

subplot(3,1,2),plot(f2,mag1,'b');  %绘出随频率变化的振幅
xlabel('频率/Hz');
ylabel('振幅');title('图2：FFT without fftshift','color','b');grid on;

subplot(3,1,3),plot(f2, mag2,'c');  %绘出随频率变化的振幅
xlabel('频率/Hz');
ylabel('振幅');title('图3：FFT after fftshift','color','c');grid on;