% Load data
function MargotData()

[~,acccalib] = load_imu('rawdata/bg14/Accelerometer/calib001.h5',[],'calib','axisorder',{'Y','Z','X'});

noise      = load_imu('rawdata/bg14/Accelerometer/noisytest001.h5');

good = noise.t >= (noise.t(end)-noise.t(1))/2;      % last half
Qgyro           = cov(noise.gyro(good,:));
constBiasGyro   = mean(noise.gyro(good,:),1);
Qacc            = cov(noise.acc(good,:));
constBiasAcc    = mean(noise.acc(good,:),1);
Qdyn            = 10*Qacc;

Qbias           = zeros(3,3); %0.05*Qgyro;
%Qbias(3,3) = 0;

imu0          = load_imu('rawdata/bg14/Accelerometer/bg14_002.h5',acccalib, ...
    'constbiasgyro',constBiasGyro, 'resamplerate',200);

knownyaw = round(linspace(1,length(imu0.t),2));
imu0.realeulerrad = zeros(size(imu0.gyro'));

imu1 = get_orient_imu(imu0, 'method','simple');
imu2 = get_orient_imu(imu0, 'method','madgwick', 'beta',0.05);
imu3 = get_orient_imu(imu0, 'method','ertss', 'Qgyro',Qgyro, ...
    'Qacc',Qacc, 'Qdyn',10*Qacc, 'Qbias',Qbias, 'Ca',1.8, 'knownorientind',knownyaw);

h(1) = subplot(4,1,1);
plot(imu1.t, imu1.orient(:,1), imu2.t,imu2.orient(:,1), imu3.t,imu3.orient(:,1));

h(2) = subplot(4,1,2);
plot(imu1.t, imu1.orient(:,2), imu2.t,imu2.orient(:,2), imu3.t,imu3.orient(:,2));

h(3) = subplot(4,1,3);
plot(imu1.t, imu1.accdyn(:,1), imu2.t,imu2.accdyn(:,1), imu3.t,imu3.accdyn(:,1));

h(4) = subplot(4,1,4);
plot(imu1.t, imu1.accdyn(:,2), imu2.t,imu2.accdyn(:,2), imu3.t,imu3.accdyn(:,2));

linkaxes(h,'x');