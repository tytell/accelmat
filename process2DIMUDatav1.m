% Vishesh Vikas
% process2DIMUDatav1(IMUfilename, encfilename, varargin)
% v1 = version 1.0
% IMUfilename = name of file containing IMU data
% encfilename = name of encoder file corresponding to IMU data
% varargin => 'static' option lets you define IMU data file for static case
% default static file = 'data_1.hdf5'
% 
% Example = process2DIMUDatav1('data_4.hdf5','enc_4','static','data_1.hdf5')

function [] = process2DIMUDatav1(IMUfilename, encfilename, varargin)
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on',...
    'DefaultAxesXMinorGrid','on','DefaultAxesYMinorGrid','on')
% Making Data Usable
%% STEP 1 - getting the biases of the sensors removed
% State of sensors - not moving, sitting on the desk
% Important = Z-Axis is pointing up

% Step 1.1 - fixing the coordinate system to body coordinate system
% Body coordinate system => y(positive) - up (gravity) and z (positive) - out of the plane
% Top accelerometer => x(positive) - up (gravity) and z (positive) - out of the plane
Rtop    = [0,-1,0;1,0,0;0,0,1];
% Bottom accelerometer => z(positive) - up (gravity) and x (positive) - right (translation)
Rbottom = [1,0,0;0,0,1;0,1,0];

% 
staticIMUfile = 'data_1.hdf5';
if ~isempty(varargin)
    for ii=1:length(varargin)
        if(strcmp(varargin{ii},'static'))
            staticIMUfile = varargin{ii+1};
        end
    end
end

% Step 2.2 - Rotating all the data to body coordinate system
staticIMUData = loadIMUData(staticIMUfile);
staticIMUData.Accel = Rtop*staticIMUData.Accel;
staticIMUData.Gyro = Rtop*staticIMUData.Gyro;
% keyboard
staticIMUData.Accel2 = Rbottom*staticIMUData.Accel2;
staticIMUData.Gyro2 = Rbottom*staticIMUData.Gyro2;

% Step 1.3 - Getting the bias of the sensors (through mean) and noise
% covariance matrix
IMU.Accelbias = mean(staticIMUData.Accel')' - [0;1;0]; % In 'g's
IMU.GyroBias = mean(staticIMUData.Gyro')';
IMU.Accelbias2 = mean(staticIMUData.Accel2')' - [0;1;0]; % In g's
IMU.GyroBias2 = mean(staticIMUData.Gyro2')';
IMU.GyroBiasComment = 'Gyroscope bias readings are in deg/sec';
IMU.AccelBiasComment = 'Accelerometer bias units are in gs not m/sec2';
IMU.AccelCovMat = cov(staticIMUData.Accel'*9.81);
IMU.Accel2CovMat = cov(staticIMUData.Accel2'*9.81);
IMU.GyroCovMat = cov(deg2rad(staticIMUData.Gyro)');
IMU.Gyro2CovMat = cov(deg2rad(staticIMUData.Gyro2)');
IMU.CovMatComment = 'Covariance matrix is in appropriate units - (rad/sec)^2 and (m/sec^2)^2';

%% STEP 2 - Remove biases and rotate readings to body coordinate system

% IMUfilename = 'data_4.hdf5';
% encfilename = 'enc_4';

theIMUData = loadIMUData(IMUfilename,'smoothwindow',200);
NN = length(theIMUData.Accel);
% keyboard
theIMUData.Accel = Rtop*theIMUData.Accel - repmat(IMU.Accelbias,1,NN);
theIMUData.Gyro = Rtop*theIMUData.Gyro - repmat(IMU.GyroBias,1,NN);
% 
theIMUData.Accel2 = Rbottom*theIMUData.Accel2 - repmat(IMU.Accelbias2,1,NN);
theIMUData.Gyro2 = Rbottom*theIMUData.Gyro2 - repmat(IMU.GyroBias2,1,NN);
figure('Name','IMU Data');
plotIMUData(theIMUData);

%% STEP 3 - Load Encoder data, smoooth it and obtain angular velocities, accelerations
theEncoder = loadEncoder(encfilename, theIMUData.deltaT,'smoothwindow',200);
plotEncoder(theEncoder);

figure('Name','Encoder vs Gyro angular velocity')
plot(theEncoder.omega);
hold on
plot(theIMUData.Gyro(3,:),'r');
title('Encoder \omega vs Gyro ');
legend('Encoder \omega','Gyroscope');
ylabel('deg/sec');

%% STEP 4 - Calcualate acceleration of the end accelerometer
d = [0;357/1000;0];  % 357mm
% theEncoder.Lambda = getLambda(deg2rad(theEncoder.omega), deg2rad(theEncoder.alpha));
theEncoder.thetaRad = deg2rad(theEncoder.theta);
theEncoder.omegaRad = [zeros(2,theEncoder.NN);deg2rad(theEncoder.omega)'];
theEncoder.alphaRad = [zeros(2,theEncoder.NN);deg2rad(theEncoder.alpha)'];
theEncoder.aBaseinBCS = zeros(3,theEncoder.NN);
for ii=1:theEncoder.NN
    ctheta = cos(theEncoder.thetaRad(ii));
    stheta = sin(theEncoder.thetaRad(ii));
    theEncoder.aBaseinBCS(:,ii) = [ctheta, stheta,0;-stheta, ctheta,0;0,0,1]*[0;9.81;0];
end
theEncoder.TopAccel = accelTwoFixedPts(theEncoder.aBaseinBCS,theEncoder.omegaRad,...
    theEncoder.alphaRad,d)/9.81;
% theEncoder.TopAccel = accelTwoFixedPts(9.81*zeros(size(theIMUData.Accel2)),theEncoder.omegaRad,...
%     theEncoder.alphaRad,d)/9.81;
figure('Name','Acceleration of top point');
subplot(2,1,1)
plot(theEncoder.TopAccel');
hold on
plot(theIMUData.Accel','--','Linewidth',2);
legend('Encoder X','Encoder Y', 'Encoder Z',...
    'IMU X', 'IMU Y', 'IMU Z');
title('Encoder vs IMU acceleration of top point');
ylabel('acceleration (g m/sec^2)');
grid on; grid minor
%% STEP 5 - Applying Extended Kalman Filter
% Xhat = EKF2Dv1(Y,u,dT,cb,ca,Qgyro,Qaccel)
% Y = accelerations
% u = gyroscope
dT = 5e-3;
cb=0.0;
ca=[0.01;0.25];
% ca=[1;0];
Xhat = EKF2Dv1(theIMUData.Accel(1:2,:)*9.81,deg2rad(theIMUData.Gyro(3,:)),dT,cb,ca,IMU.GyroCovMat,IMU.AccelCovMat);
% figure('Name','Encoder vs Estimated')
subplot(2,1,2)
plot(rad2deg(Xhat(1,:)));
hold on; plot(theEncoder.theta)
legend('Estiamted','Encoder');
ylabel('\theta (deg)');
grid on; grid minor;
title('\theta Encoder vs Estimated');
end

%%
function imuData = loadIMUData(filename, varargin)
    imuData = h52Struct(filename);
    % 
    imuData.smoothwindow = 1;
    if ~isempty(varargin)
        for ii=1:2:length(varargin)
            if(strcmp(varargin{ii},'smoothwindow'))
                imuData.smoothwindow = varargin{ii+1};
            end
        end
    end
    % Step 1 - Gathering data
    % Time stamp on the data
    imuData.timestamp = imuData.t';
    % The frequency of observation
    imuData.deltaT = 5e-3;
    % The actual time in milliseconds
    imuData.t = imuData.timestamp*imuData.deltaT;
    imuData.comment = 'suffix 2 corresponds to base IMU';
    imuData.AccelNorm = vectorNorm(imuData.Accel);
    imuData.Accel2Norm = vectorNorm(imuData.Accel2);
    % Step 2 - smooth data
    for ii=1:3
        imuData.Accel(ii,:)   = smooth(imuData.Accel(ii,:),imuData.smoothwindow);
        imuData.Gyro(ii,:)    = smooth(imuData.Gyro(ii,:),imuData.smoothwindow);
        imuData.Accel2(ii,:)  = smooth(imuData.Accel2(ii,:),imuData.smoothwindow);
        imuData.Gyro2(ii,:)   = smooth(imuData.Gyro2(ii,:),imuData.smoothwindow);
    end
end

%%
% Vishesh Vikas
% function encoder = loadEncoder(filename,deltaT,varargin)
function encoder = loadEncoder(filename,deltaT,varargin)
    encoder.smoothwindow    = 100;
    for ii=1:2:length(varargin)
        if strcmp(varargin{ii},'smoothwindow')
            encoder.smoothwindow    = varargin{ii+1};
        end
    end
    aa = load(filename);
    encoder.t               = (0:(length(aa)-1))'*deltaT;
    encoder.NN              = length(aa);
    encoder.theta           = smooth(aa,encoder.smoothwindow);
    encoder.theta           = encoder.theta - ones(size(encoder.theta))*encoder.theta(1);
    encoder.omega           = smooth(diff(aa),encoder.smoothwindow)/deltaT;
    encoder.alpha           = smooth(diff(encoder.omega),encoder.smoothwindow)/deltaT;
    encoder.omega           = [encoder.omega(1);encoder.omega];
    encoder.alpha           = [encoder.alpha(1);encoder.alpha;encoder.alpha(end)];
end

%%
function plotIMUData(imuData)
    % Step 2 - Plotting datafigure()
    subplot(2,2,1)
    plot(imuData.t,imuData.Accel'); 
    title('Accelerometer 1 (Top)');
    xlabel('time (ms)');
    ylabel('g (9.81 m/s^2)');
    legend('x','y','z');
    subplot(2,2,2)
    plot(imuData.t,imuData.Accel2')
    title('Accelerometer 2 (Base)');
    xlabel('time (ms)');
    ylabel('g (9.81 m/s^2)');
    legend('x','y','z');    
    subplot(2,2,3)
    plot(imuData.t,imuData.Gyro'); 
    title('Gyroscope 1 (Top)');
    xlabel('time (ms)');
    ylabel('deg/sec');
    legend('x','y','z');
    subplot(2,2,4)
    plot(imuData.t,imuData.Gyro2')
    title('Gyroscope 2 (Base)');
    xlabel('time (ms)');
    ylabel('deg/sec');
    legend('x','y','z');
end

%%
% Vishesh Vikas
% function encoder = plotEncoder(encoder)
% encoder = struct obtained from loadEncoder(filename, deltaT)
function plotEncoder(encoder)
figure();
subplot(3,1,1)
% keyboard
plot(encoder.t,encoder.theta);
% xlabel('time (ms)');
ylabel('angle \theta (deg)');
xlabel('time (ms)');
grid on;
subplot(3,1,2)
plot(encoder.omega);
ylabel('anglular velocity \omega (deg/sec)');
xlabel('time (ms)');
grid on;
subplot(3,1,3)
plot(encoder.alpha);
ylabel('anglular acceleration \alpha (deg/sec^2)');
xlabel('time (ms)');
grid on;
end

%%
% Vishesh Vikas
% aA,omega,alpha = 3XN, rAB = 3X1
function aB = accelTwoFixedPts(aA,omega,alpha,rAB)
    crossMat = @(x) [0,-x(3),x(2); x(3),0,-x(1); -x(2),x(1),0];
%     lambda = getLambda(omega,alpha);
    aB = zeros(size(aA));
    for ii=1:length(omega)
        aB(:,ii) = aA(:,ii) + (crossMat(alpha(:,ii)) + ...
            crossMat(omega(:,ii))^2)*rAB;
    end
end

%%
% Vishesh Vikas
% function Xhat = EKF2Dv1(Y,u,dT,cb,ca,Qgyro,Qaccel)

function Xhat = EKF2Dv1(Y,u,dT,cb,ca,Qgyro,Qaccel)
% System Dynamics
% X(k) = [theta;bGyro;aDyninN]
% theta_Dot = -bGyro + omegaGyro - eGyro
% bGyro_Dot = cb
% aDyninN_Dot = ca

% Discrete-time Algebraic Reccati Equation based EKF
% X(k+1) = f(x)
% y(k) = h(x)
% ca = ca*ones(2,1);
NN      = length(Y);
QX      = 0.1*eye(4);
% QX = rand(4); QX=0.1*QX'*QX;
Qk      = diag([Qgyro(1,1), zeros(1,3)])*dT^2 + QX;
Rk      = Qaccel(1:2,1:2);
% SkData      = zeros(2,2,NN+1); % Size of Rk
Pkgkm1Data  = zeros(4,4,NN);
% PkgkData    = zeros(4,4,NN);
xhatkgkm1Data = zeros(4,NN+1);
% xhatkgkData = zeros(4,NN);
Pkgkm1Data(:,:,1) = eye(4);
xhatkgkm1Data(:,1) = zeros(4,1);
for ii=1:NN
    Xk  = xhatkgkm1Data(:,ii);
    uk  = u(:,ii);
    yk  = Y(:,ii);
%     [fk,fkPrime,fkDPrime] = getF(Xk,dT,uk,cb,ca);
    [hk, hkPrime]   = getH(Xk);
    Pkgkm1          = Pkgkm1Data(:,:,ii);
    xhatkgkm1       = xhatkgkm1Data(:,ii);
%     
    Sk              = Rk + hkPrime*Pkgkm1*transpose(hkPrime);
    Kk              = Pkgkm1*transpose(hkPrime)*pinv(Sk);
    ek              = yk - hk;
    xhatkgk         = xhatkgkm1 + Kk*ek;
    Pkgk            = Pkgkm1 - Kk*hkPrime*Pkgkm1;
    [fk,fkPrime,fkDPrime] = getF(xhatkgk,dT,uk,cb,ca);
    xhatkgkm1Data(:,ii+1)   = fk;
    Pkgkm1Data(:,:,ii+1)    = Qk + fkPrime*Pkgk*transpose(fkPrime);
end
Xhat = xhatkgkm1Data(:,2:end);
end

function [fk,fkPrime,fkDPrime] = getF(Xk,dT,uk,cb,ca)
    JX  = eye(4);
    JX(1,2) = -dT;
    fk  = JX*Xk+ ...
        [0;cb;ca]*dT +[dT;0;0;0]*uk;
    fkPrime = JX;
    fkDPrime = zeros(3,3);
end

function [hk, hkPrime] = getH(Xk)
    theta   = Xk(1);
    ct      = cos(theta);
    st      = sin(theta);
    Rk      = [ct, st;-st,ct];
    RkPrime = [-st, ct;-ct,-st];
    aDyn    = Xk(3:4);
    gN      = [0;1];
    hk      = Rk*(gN +aDyn);
    hkPrime = [RkPrime*(gN +aDyn), zeros(2,1), Rk];
end

%%
% Vishesh Vikas - vishesh.vikas@gmail.com
% 2nd norm of vector 'vec'
% input = vec 3 X N 
% output = normV 1 X N
% normV(i) = ||vec(:,i)||_2 (2nd norm)
function normV = vectorNorm(vec)
normV = sqrt(sum(vec.*vec,1));
end