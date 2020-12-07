function [scanpattern, doas] = MyMusic2D(X, pos, f, c, Aziang, eleang, NumSignals)

if nargin<7
    % 用AIC或者MDL的方法确定信号数量
    NumSignals = 1;
end

% 计算协方差矩阵
T = size(X);
Cx = X'*X/T(1);
% 计算特征值和特征向量，并从小到大排序
[V, D] = eig(Cx);
[eigenvals,i] = sort(diag(D));
eigenvects = V(:,i);
% 得到噪声对于的特征向量
noise_eigenvects = eigenvects(:,1:end-NumSignals);
% 计算迭代次数，最小化计算开销
len_az = numel(Aziang);
len_el = numel(eleang);
patternSize = [len_el, len_az];
[scanAz, scanEl] = meshgrid(Aziang,eleang);
scanAngles = [scanAz(:) scanEl(:)].';
numIter = min(len_az, len_el); % 循环次数
scanAngleBlockSize = max(len_az, len_el); % 每次循环扫描多少个角度
scanAngleBlockIndex = 1:scanAngleBlockSize;
% scan
for idx = 0:numIter-1
    curIdx = scanAngleBlockIndex + idx*scanAngleBlockSize;
    curAngle = scanAngles(:, curIdx);
    % 计算steer矩阵A
    a = steer_matrix(pos, f, c, curAngle);
    % 等同于
    % A = sv'*noise_eigenvects*noise_eigenvects'*sv
    % D = diag(A,0) + eps(1)
    % 这样写减少计算开销
    D = sum(abs((a'*noise_eigenvects)).^2,2)+eps(1);
    pattern(curIdx) = 1./D;
end
classtouse = class(X);
scanpattern = cast(reshape(sqrt(pattern),patternSize),classtouse);
doas = 0;
end

function a = steer_matrix(pos, f, c, scanAngles)
% az是与x轴的夹角，el是与xoy平面的夹角
azang = scanAngles(1,:);
elang = scanAngles(2,:);
% 前提应该是第一个麦克风坐标在圆心才行，为什么matlab不是?
% 猜测在仿真时，sensorsig()函数是对着原点进行延迟构造的，所以仿真没问题，放到真实环境有问题。
signalDir = [-cosd(elang).*cosd(azang);...
    -cosd(elang).*sind(azang);...
    -sind(elang)];
timeDelay = pos.'*signalDir/c;
a = exp(-1i*2*pi*f*timeDelay);
end