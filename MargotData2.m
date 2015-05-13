% Load data
function MargotData2()

theNoise01      = load_imu('rawdata/bg14/Accelerometer/noisytest001.h5');
theNoise01.gyro = deg2rad(theNoise01.gyro);
constBiasGyro   = mean(theNoise01.gyro(1.5e4:8e4,:),1);
constBiasAcc    = mean(theNoise01.acc(100:8e4,:),1);

Qgyro           = cov(theNoise01.gyro(1.5e4:8e4,:));
Qacc            = cov(theNoise01.acc(100:8e4,:));
QaccAmp         = 1;
Qdyn            = QaccAmp*Qacc;
Qbias           = diag([0.01,0.01,1e-5])*Qacc;
Ca              = 0.9;

[~,acccalib] = load_imu('rawdata/bg14/Accelerometer/calib001.h5',[],'calib','axisorder',{'Y','-Z','-X'});

theData         =  load_imu('rawdata/bg14/Accelerometer/bg14_004.h5');

theData.gyro    = deg2rad(theData.gyro);

NN              = length(theData.t);
zeroInd         = find(theData.t==0);
zeroNN          = length(zeroInd);


sN              = zeroInd(1)+1;
endN            = zeroInd(2)-1;
theData.gyro    = theData.gyro(sN:endN,:);
theData.acc     = theData.acc(sN:endN,:);
theData.t       = theData.t(sN:endN);
NN              = length(theData.t);

theData.gyro    = theData.gyro - repmat(constBiasGyro,NN,1);
% 
knownyaw = [1 length(theData.gyro)];
theData.realeulerrad = zeros(size(theData.gyro'));

results         = ERTSSv1(theData, Qgyro, Qbias, Qacc, Qdyn, Ca, knownyaw);

% figure('Name','The data');
figure('name', strcat('Ca=',num2str(Ca),' QaccAmp=',num2str(QaccAmp)));
subplot(2,2,1)
plot(theData.t,theData.acc);
title('Data - acceleration');
legend('X','Y','Z');
xlabel('time'); ylabel('acc (m/s^2)');
subplot(2,2,3)
plot(theData.t, rad2deg(theData.gyro));
title('Data - gyroscope');
legend('X','Y','Z');
xlabel('time'); ylabel('angular vel (deg/s)');
% 
subplot(2,2,2);
% plot(rad2deg(eulerEFK'),'-.','Linewidth',2);
% hold on
plot(theData.t,rad2deg(results.euler'),'Linewidth',2);
title('Filtered Orientation')
legend('Roll','Pitch','Yaw');
xlabel('time'); ylabel('angle (degrees)');
subplot(2,2,4);
% keyboard
plot(theData.t,results.dynamicAcc')
title('Dynamic Acceleration');
legend('X','Y','Z');
xlabel('time'); ylabel('acc (m/s^2)');


end