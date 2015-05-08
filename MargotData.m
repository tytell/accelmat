% Load data
function MargotData()

theNoise01      = load_imu('rawdata/bg14/Accelerometer/noisytest001.h5');
theNoise01.gyro = deg2rad(theNoise01.gyro);
Qgyro           = cov(theNoise01.gyro(1.5e4:8e4,:));
constBiasGyro   = mean(theNoise01.gyro(1.5e4:8e4,:),1);
Qacc            = cov(theNoise01.acc(100:8e4,:));
constBiasAcc    = mean(theNoise01.acc(100:8e4,:),1);
Qdyn            = 10*Qacc;

Qbias           = zeros(3,3);

[~,acccalib] = load_imu('rawdata/bg14/Accelerometer/calib001.h5',[],'calib','axisorder',{'Y','-Z','-X'});

theData         =  load_imu('rawdata/bg14/Accelerometer/bg14_004.h5');
%theData         = load('Bg5-20_0hz-005.mat');
theData.gyro = deg2rad(theData.gyro);
% keyboard
% [results]       = smootherRTS2(theData, Qgyro, Qbias, Qacc);
% results         = smootherRTSwithAv3(theData, Qgyro, Qbias, Qacc, 10*Qdyn, -0.75, 10);

knownyaw = [1 length(theData.gyro)];
theData.realeulerrad = zeros(size(theData.gyro'));

results         = ERTSSv1(theData, Qgyro, Qbias, Qacc, Qdyn, 1.1, knownyaw);
