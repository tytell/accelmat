% Load data
function MargotData()

theNoise01      = load_imu('rawdata/bg14/Accelerometer/noisytest001.h5');
theNoise01.gyro = deg2rad(theNoise01.gyro);
Qgyro           = cov(theNoise01.gyro(1.5e4:8e4,:));
constBiasGyro   = mean(theNoise01.gyro(1.5e4:8e4,:),1);
Qacc            = cov(theNoise01.acc(100:8e4,:));
constBiasAcc    = mean(theNoise01.acc(100:8e4,:),1);
Qdyn            = 10*Qacc;

Qbias           = 0.05*Qgyro;

[~,acccalib] = load_imu('rawdata/bg14/Accelerometer/calib001.h5',[],'calib','axisorder',{'Y','-Z','-X'});

theData         =  load_imu('rawdata/bg14/Accelerometer/bg14_015.h5');
% Converting from degrees to radians for calculations
theData.gyro = deg2rad(theData.gyro);
% 

theDattheData.gyro    = theData.gyro - repmat(constBiasGyro,NN,1);

% subtracting the constant bias from gyro readings
theData.gyro    = theData.gyro - repmat(constBiasGyro,NN,1);
good = theData.t >= -20;
theData.gyro = theData.gyro(good,:);
theData.acc = theData.acc(good,:);
theData.t = theData.t(good,:);

knownyaw = round(linspace(1,length(theData.t),2));
theData.realeulerrad = zeros(size(theData.gyro'));

results         = ERTSSv1(theData, Qgyro, Qbias, Qacc, Qdyn, 1.1, knownyaw);
