% This function generates simulated IMU data of a virtual fish

function [simimu] = simluatedData()

samplefreq	= 1e-3;     % 1 KHz
time        = 0:samplefreq:10; % in seconds
% Accelerometer (datesheet) quantities
accel.noisestd  = (3e-4*sqrt(0.5/samplefreq));          % unit = g (9.81 m/s2)
% Gyroscope (datasheet) quantities
gyro.noisestd   = deg2rad(0.01*sqrt(0.5/samplefreq)); % unit = radians
gyro.biasstd    = 2*gyro.noisestd;
% Simulation characteristics
freq        = 0.5;

roll.amplitude  = deg2rad(10);
roll.phase      = deg2rad(25);
roll.freq       = freq;
roll.angle      = roll.amplitude*sin(2*pi*roll.freq*time + roll.phase);
roll.rate       = 2*pi*roll.freq*roll.amplitude*cos(2*pi*roll.freq*time + roll.phase);
roll.rateofrate = -(2*pi*roll.freq)^2*roll.angle;

pitch.amplitude = 0;
pitch.angle     = zeros(size(time));
pitch.rate      = zeros(size(time));
pitch.rateofrate = zeros(size(time));

yaw.amplitude   = deg2rad(20);
yaw.phase       = deg2rad(40);
yaw.freq        = freq;
yaw.angle       = yaw.amplitude*sin(2*pi*yaw.freq*time + yaw.phase);
yaw.rate        = 2*pi*yaw.freq*yaw.amplitude*cos(2*pi*yaw.freq*time + yaw.phase);
yaw.rateofrate  = -(2*pi*yaw.freq)^2*yaw.angle;

% Converting Roll-Pitch-Yaw into Euler angles
% RPY = ZYX Euler Angles = R(x,phi)*R(y,theta)*R(z,psi)
% psi = Yaw, theta = Pitch, phi = Roll
psi     = yaw;
theta   = pitch;
phi     = roll;
 [Omega_Sensor, Alpha_Sensor] = ...
    getRotationRatesEulerZYX(psi, theta, phi);


% Rotation of the rigid body
Ax          = 3*(accel.noisestd);
rOP         = [0,0,0]'; % meters
Accel_CoM     = [Ax*sin(4*pi*freq*time);zeros(2,length(time))];
Accel_P     = getAcceleration(Accel_CoM, rOP, Omega_Sensor, Alpha_Sensor);

% Adding noise
% gyro.sigma = (gyro.noisevar)*eye(3,3); R = (chol(gyro.sigma));
simimu.gyro     = Omega_Sensor + (gyro.noisestd)*randn(size(Omega_Sensor));
% Adding bias drift
Bias            = getBias(length(Omega_Sensor), gyro.biasstd);
simimu.gyro     = simimu.gyro+Bias;
% accel.sigma = (accel.noisevar)*eye(3,3); R = chol(accel.sigma);
simimu.acc  = Accel_P + accel.noisestd*randn(size(Accel_P));

simimu.gyro     = simimu.gyro';
simimu.acc      = simimu.acc';
simimu.t        = time';


figure('Name','Gyroscope with noise, without drift')
plot(time, rad2deg(Omega_Sensor));
hold on
plot(time, rad2deg(simimu.gyro));
legend('Sensor X', 'Sensor Y', 'Sensor Z')

figure()
plot(time, Accel_P);
hold on
plot(time, simimu.acc);

end

function [Omega_Sensor, Alpha_Sensor] = ...
    getRotationRatesEulerZYX(psi, theta, phi)
% 
% keyboard
    N = length(psi.angle);
    Omega_Sensor = zeros(3,N);
    Alpha_Sensor = zeros(3,N);
    for ii=1:N
        Theta_d     = [phi.rate(ii); theta.rate(ii); psi.rate(ii)];
        Theta_dd    = [phi.rateofrate(ii); theta.rateofrate(ii); psi.rateofrate(ii)];
        W = [1, 0, -sin(theta.angle(ii));
            0, cos(phi.angle(ii)), cos(theta.angle(ii)).*sin(phi.angle(ii));
            0, -sin(phi.angle(ii)), cos(theta.angle(ii)).*cos(phi.angle(ii))];
        W_theta =[0,0,-cos(theta.angle(ii));
            0,0,-sin(theta.angle(ii))*sin(phi.angle(ii));
            0,0,-sin(theta.angle(ii))*cos(phi.angle(ii))];
        W_phi   = [0,0,0;
            0,-sin(phi.angle(ii)), cos(theta.angle(ii))*cos(phi.angle(ii));
            0,-cos(phi.angle(ii)), -cos(theta.angle(ii))*sin(phi.angle(ii))];
        W_d     = W_theta*theta.rate(ii) + W_phi*phi.rate(ii);
        Omega_Sensor(:,ii)     = W*Theta_d;
        Alpha_Sensor(:,ii)     = W_d*Theta_d + W*Theta_dd;
    end
end

% Two points on a rigid body
function [Accel_P] = getAcceleration(Accel_O, rOP, Omega_Sensor, Alpha_Sensor)
    N = length(Accel_O);
    Accel_P = zeros(3,N);

    for ii=1:N
        Accel_P(:,ii) = Accel_O(:,ii) + ...
            (crossMat(Alpha_Sensor(:,ii)) + ...
            crossMat(Omega_Sensor(:,ii))^2)*rOP;
    end
end

% Cross Product matrix
function [Mat] = crossMat(vector)
    Mat = [0, -vector(3), vector(2);
        vector(3), 0, -vector(1);
        -vector(2), vector(1), 0];
end

% Random walk
function [Bias] = getBias(length, std)
    Bias = zeros(3,length);
    Bias(:,1) = std*randn(3,1);
    for ii=2:length
        Bias(:,ii) = Bias(:,ii-1) + std*randn(3,1);
    end
end