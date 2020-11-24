dt=0.00001; %设定步长

t=0:dt:3; %设置3个频率的正弦信号 300HZ ,500HZ,1000HZ

s1=sin(2*pi*300*t);

s2= sin(2*pi*500*t);

s3= sin(2*pi*1000*t);

s=s1+s2+s3; % 3个正弦信号叠加

figure(1);

subplot(2,2,1);plot(t,s1);

xlabel('t');title('300HZ正弦信号');

set(gca, 'XLim', [0,0.01]);

subplot(2,2,2); plot(t,s2);

xlabel('t');title('500HZ正弦信号');
set(gca, 'XLim', [0,0.01]);

subplot(2,2,3); plot(t,s3);

xlabel('t');title('1000HZ正弦信号');
set(gca, 'XLim', [0,0.01]);

subplot(2,2,4); plot(t,s);

xlabel('t');title('合成信号');
set(gca, 'XLim', [0,0.01]);


ss=fft(s,4096);

SS=(abs(ss(1:1:2049))); %求合成信号频谱

k1=0:2048;

w1=(1/4096)*k1*10000; %取0......Fs/2的部分

figure(2);

plot(w1,SS); grid %画频谱图

title('求原信号频谱');

%****通过低通滤波器*****%

ws1=1000; %设计一个带通为300HZ，阻带为1000HZ的低通滤波器

wp1=300; wc=5000;

wp=wp1/wc; ws=ws1/wc;

[n,wn]=buttord(wp,ws,1,30) %巴特沃斯低通滤波器

[b,a]=butter(n,wn);

sb=3*filter(b,a,s) ; %合成信号通过低通滤波器

ssb=fft(sb,4096); %求频谱

SSb=abs(ssb(1:1:2049));

k1=0:2048; w1=(1/4096)*k1*10000; %画频谱图

figure(3);

plot(w1,SSb); grid

title('经过低通滤波器后的信号频谱');