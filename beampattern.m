clc;clear;
close all;
numElements = 10;      
spacing = 0.2;      
freq = 1000;
speedSound = 343.0;   
ANGLE_RESOLUTION = 500;

r = zeros(500,1);
theta = zeros(500,1);

% x取什么值都一样，很神奇
x = 0;

for a=0:1:ANGLE_RESOLUTION-1
    angle = 0 + 360.0 * a / (ANGLE_RESOLUTION-1);
    angleRad = pi * angle / 180;
    realSum = 0;
    imagSum = 0;
    for i=0:1:numElements-1
        position = i * spacing;
        delay = position * sin(angleRad) / speedSound;
%         if i == 0
%             delay = 0;
%         else
%             delay = 0.043 * cos((i-1)*(1/3*pi)-angleRad);
%         end
        realSum = realSum + cos(2.0 * pi * freq * delay + x);
        imagSum = imagSum + sin(2.0 * pi * freq * delay + x);
    end
    
    output = sqrt(realSum * realSum + imagSum * imagSum) / numElements;
    logOutput = 20 * log10(output);
    if (logOutput < -50) 
        logOutput = -50;
    end
    r(a+1) = logOutput;
    theta(a+1) = angleRad;
end
polarplot(theta,r);
rlim([-50 0])