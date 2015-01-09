% Kamlman Filter for measuring orientation using gyroscope and
% accelerometer
% 
% Research paper to follow - 
% "Measuring orientation of human body segments using miniature guroscopes
% and acclerometers", H.J. Luinge, P.H. Veltink, Medical & Biological
% Engineering & Computin
% 
% Author - Vishesh Vikas (vishesh.vikas@gmail.com)
% 
% IMU contains one tri-axial acceleromter and one triaxial gyroscope
% Sensor Fusion with a Kalman Filter -
% Zh_G and Zh_A are estimates of gravitional acceleration/inclination from
% gyroscope and acclerometer but with different accuracy
function [Xk] = KalmanVeltink(imu)
% imu.gyro = N x 3, imu.accel = N x 3

systemDynamics.A = zeros(6,6);
systemDynamics.B = 0;


% bt0 = mean(imu.gyro(1:2000,:))';
bt0     = imu.gyro(1,:)';
aGt0    = imu.acc(1,:)';
R0      = eye(3);
N       = size(imu.gyro,1);
g       = 9.81;
Ca      = 0.3;
QvG     = 0.2*eye(3); % gyroscope reading
QvA     = 0.2*eye(3); % accelerometer reading
QwA     = 0.01*eye(3); % acceleration bias drift 
Qb      = 0.01*eye(3);
Qtheta_tm1p = eye(3);
Qb_tm1p = eye(3); 
Qa_tm1p = eye(3);
b_tm1p  = bt0;
aG_tm1p = aGt0;
RotMat_tm1p = R0;

Pkm1 = eye(6);
xkm1 = zeros(6,1);
% ====> Initializing 
% theErr = 100;
    for ii=2:N
%     while(abs(theErr)>1e-2)
%         systemDynamics.A = zeros(6,6);
%         systemDynamics.B = 0;
    yGt = imu.gyro(ii,:)';
    yAt = imu.acc(ii,:)'*g;
    T   = imu.t(ii)-imu.t(ii-1);
    [ZGyro_tm, ZAccel_tm, aS_tm] = incliantionEstimates(yGt, yAt, aG_tm1p, ...
        RotMat_tm1p, b_tm1p, T, Ca);
    systemDynamics.C = errorDynamics(ZAccel_tm, aS_tm, ZGyro_tm, g, T);
    [Qwt, Qvt] = covarianceMatricies(Qtheta_tm1p, Qb_tm1p, Qa_tm1p, ...
    QvG, QvA, QwA, Qb, g ,T, Ca);
    systemDynamics.Q = Qwt;
    systemDynamics.R = Qvt;
    zk = ZAccel_tm - ZGyro_tm;
    uk = 0;
    [xk, Pk] = theKalmanFilterStep(systemDynamics, uk, Pkm1, xkm1, zk);
    Xk(:,ii) = xk;
%     resetting old values
    Pkm1 = Pk;
    Qtheta_tm1p = Pk(1:3,1:3);
    Qb_tm1p     = Pk(4:6,4:6);
    xkm1 = xk;
    end
end


% ---------------------------------------------------------------------
% STEP 1: Obtain inclination/gravity in sensor coordinate system (S) from
% gyroscope and accelerometer

% Inputs required are
% Gyro Signal           - yGt
% Accelrometer Signal   - yAt
% Linear acceleration estimate in Global CS for t-1 (tm1) - aG_tm1p
% Rotation matrix estimate for t-1 (tm1) - RotMat_tm1p
% Gyroscope bias estiamte for t-1 (tm1) - b_tm1p
% Time step             - T
% Acceleration change (low pass filter) constant - Ca
% Non-gravity acceleration - aS_tm
function [ZGyro_tm, ZAccel_tm, aS_tm] = incliantionEstimates(yGt, yAt, aG_tm1p, ...
    RotMat_tm1p, b_tm1p, T, Ca)
% Step 1a: Get angular velocity estimate omega_tm
omega_tm = yGt - b_tm1p;
% Step 1b: Get estimate for Rotation matrix
RotMat_tm = RotMat_tm1p + RotMat_tm1p*crossProductMatrix(omega_tm)*T;
% Step 1c: Get gravity vector estimate in S using gyroscope
ZGyro_tm = RotMat_tm(3,:)'; % Column vector
% Step 1d: non-gravity linear acceleration estimate
aG_tm = Ca*aG_tm1p;
% Step 1e: non-gravity linear acceleration in S
aS_tm = RotMat_tm*aG_tm;
% Step 1f: Gravity vector estimate in S using Accelerometer
ZAccel = yAt - aS_tm;
ZAccel_tm = ZAccel/norm(ZAccel);
end


% Matrix Equivalent representation of cross product
% theMatrix = Equivalent cross product matrix
% theVector = cross product vector
% theMatrix is skew symmetric
function theMatrix = crossProductMatrix(theVector)
theMatrix = [0 -theVector(3), theVector(2);
    theVector(3), 0, theVector(1);
    -theVector(2), -theVector(1), 0];
end

% The error dynamics
% x(k+1) = A*x(k) + wt
% y(k+1) = C*x(k) + vt
% x(k) = [theta_error; bias_error]
% theta_error = 3x1 because it's the change in axis of rotation
% Z_tm1 = Z estimate at t-1
% aS_tm = non-gravity acceleration at t before filtering
% T = time change (delta T)
function [C] = errorDynamics(ZA_tm1, aS_tm, ZG_tm1, g, T)
% A = 0;
C = [crossProductMatrix(ZA_tm1- aS_tm/g),crossProductMatrix(T*ZG_tm1)];
end
% function [C] = errorDynamics(Z_tm1, aS_tm, g, T)
% % A = 0;
% C = crossProductMatrix(Z_tm1- aS_tm/g)*crossProductMatrix(T*Z_tm1);
% end

% Covariance matricies
% QvG = Gyroscope noise covariance
% QvA = Accelerometer noise covariance
% QwA = Acceleration rate of change covariance (low pass filter)
% Qb = Gyroscope bias noise covariance
% Equation 24 in the paper
function [Qwt, Qvt] = covarianceMatricies(Qtheta_tm1p, Qb_tm1p, Qa_tm1p, ...
    QvG, QvA, QwA, Qb, g ,T, ca)
Qwt = [Qtheta_tm1p+ T^2*Qb_tm1p + T^2 QvG;
    T^2*Qb_tm1p, Qb_tm1p + Qb];
Qvt = QvG + (ca^2*Qa_tm1p + QwA + QvA)/g^2;
end

% The Kalman Filter step
% 
% The error dynamics - z_error = ZAccel - ZGyro
% x(k+1) = A*x(k) + B*u(k) + w(k)
% z(k+1) = C*x(k) + v(k)
% w = N(0,P)
% v = N(0,Q)
% Figure 4.2 from Welch-Bishop Kalman paper
function [xk, Pk] = theKalmanFilterStep(systemDynamics, uk, P_km1, xkm1, zk)
% Step 1 : Time Update
% x_prior(k) = A*x_posterior(k-1) + B*u(k)
xk_prior = systemDynamics.A*xkm1 +  systemDynamics.B*uk;
% Error dynamics covariance
%P_prior(k) = A*P_posterior(k-1)*A' + Q
Pk_prior = systemDynamics.A*P_km1*systemDynamics.A' + systemDynamics.Q;
% Step 2 : Measurement Update
% Kalman Gain
% keyboard
Kk = Pk_prior*systemDynamics.C'*...
    inv(systemDynamics.C*Pk_prior*systemDynamics.C' + systemDynamics.R);
% measurement update
xk = xk_prior + Kk*(zk-systemDynamics.C*xk_prior);
% covariance matrix update
Pk = (eye(6)-Kk*systemDynamics.C)*Pk_prior;
end
