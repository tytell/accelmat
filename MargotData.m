% Load data
function MargotData()

[cal,acccalib] = load_imu('rawdata/bg14/Accelerometer/calib001.h5',[],'calib','axisorder',{'-Y','-Z','X'});

noise      = load_imu('rawdata/bg14/Accelerometer/noisytest001.h5');

good = noise.t >= (noise.t(end)-noise.t(1))/2;      % last half
Qgyro           = cov(noise.gyro(good,:));
constBiasGyro   = mean(noise.gyro(good,:),1);
Qacc            = cov(noise.acc(good,:));
constBiasAcc    = mean(noise.acc(good,:),1);
Qdyn            = 10*Qacc;

Qbias           = zeros(3,3); %0.05*Qgyro;
%Qbias(3,3) = 0;

beta = sqrt(3/4) * 2*max(rms(noise.gyro(good,:)));

imu0          = load_imu('rawdata/bg14/Accelerometer/bg14_004.h5',acccalib, ...
    'constbiasgyro',constBiasGyro, 'resamplerate',200);

knownyaw = round(linspace(1,length(imu0.t),2));
imu0.realeulerrad = zeros(size(imu0.gyro'));

imu1 = get_orient_imu(imu0, 'method','simple');
imu2 = get_orient_imu(imu0, 'method','madgwick', 'beta',beta);
%imu3 = get_orient_imu(imu0, 'method','ertss', 'Qgyro',Qgyro, ...
%    'Qacc',Qacc, 'Qdyn',10*Qacc, 'Qbias',Qbias, 'Ca',1.0, 'knownorientind',knownyaw);

h(1) = subplot(5,1,1);
plot(imu1.t, imu1.orient(:,1), imu2.t,imu2.orient(:,1)); %, imu3.t,imu3.orient(:,1));
ylabel('Roll');
legend('simple','madgwick','ertss');

h(2) = subplot(5,1,2);
plot(imu1.t, imu1.orient(:,2), imu2.t,imu2.orient(:,2)); %, imu3.t,imu3.orient(:,2));
ylabel('Pitch');

h(3) = subplot(5,1,3);
plot(imu1.t, unwrap(imu1.orient(:,3)), imu2.t,unwrap(imu2.orient(:,3))); %, imu3.t,imu3.orient(:,3));
ylabel('Yaw');

h(4) = subplot(5,1,4);
plot(imu1.t, imu1.accdyn(:,1), imu2.t,imu2.accdyn(:,1)); %, imu3.t,imu3.accdyn(:,1));
ylabel('Forward acc');

h(5) = subplot(5,1,5);
plot(imu1.t, imu1.accdyn(:,2), imu2.t,imu2.accdyn(:,2)); %, imu3.t,imu3.accdyn(:,2));
ylabel('Side acc');

linkaxes(h,'x');