%% ERTSS for simulation data
% [results] = simERTSS(imu, Ca, QaccAmp)
% imu = contains gyro signal, accelerometer signal, etc.
% Ca = constant for low pass filter of dynamic acceleration.
function [results] = simERTSS3(imu, Ca, QaccAmp)
% Step 1 - get all data
Gyro        = imu.gyro'; 
Acc         = imu.acc';
dT          = diff(imu.t);
dT          = [dT(1); dT];
% Step 2a - setting up all the covariance matrices
Qacc        = QaccAmp*imu.accnoisestd*eye(3);
Qdyn        = imu.accnoisestd*eye(3);
Qgyro       = imu.gyronoisestd*eye(3);
Qbias       = 0.1*imu.gyrobiasdriftstd*eye(3);

% Observation covariance matrix Rk - constant as we observe the
% accelerometer
Rk      = Qacc;

% Initializing the Extended Rauch-Tung-Strielbel smoother (ERTSS)
xkm1      = zeros(9,1);
T         = dT(1);
Pkm1      = [(Qgyro+Qbias)*T^2, -Qbias*T, zeros(3,3);-Qbias*T, Qbias,zeros(3,3);...
    zeros(3,6), Qdyn];
% Gravity - magnitude is 1 and sign is negative
gN      = [0,0,-1]';
gN = gN/norm(gN);

% Defining matrices that store calculated data
NN          = length(Gyro);
eulerEFK    = zeros(3,NN);      % euler angles post filtering
eulerS      = zeros(3,NN);      % euler angles post smoothing
PkEKF       = zeros(9,9,NN);    % innovation covariance matrices
err         = zeros(3,NN);      %

% the data point where 'yaw' is known
knownN = sort(floor(random('unif',1,NN,20,1)));
knownN = [1;knownN;NN];
% knownN = [1;NN];

tic
% Step 3 - FORWARD FILTERING STEP
for ttfwd=1:NN
    T = dT(ttfwd);
%     Gyrosocpe reading and accelerometer readings
    Omegak          = Gyro(:,ttfwd);
    Accel           = Acc(:,ttfwd);
%     Step 3a - system dynamics and conditional gaussian
    [Fk,xkM,Bk]     = systemDynamics(xkm1, Omegak, T,Ca);
    theQT           = getQT(xkM(1:3));
    Qk              = [Bk*(Qgyro+Qbias)*Bk'*T^2, -Bk*Qbias*T, zeros(3,3);...
        -Qbias*Bk'*T, Qbias, zeros(3,3);
        zeros(3,6), Qdyn];
    PkM             = Fk*Pkm1*Fk' + Qk;
    if(sum(ttfwd==knownN)==1)
        [hk2, Jh2] = observationDynamicsEKF2(xkM, gN);
        Hk          = Jh2;
        Sk              = Hk*PkM*Hk' + [Rk,zeros(3,1);zeros(1,4)];
        Kk              = PkM*Hk'*pinv(Sk);
        xk          = xkM + Kk*([Accel;imu.realeulerrad(3,ttfwd)]...
                                -hk2);
    else 
        [hk, Jh] = observationDynamicsEKF(xkM, gN);
        Hk = Jh;
        Sk              = Hk*PkM*Hk' + Rk;
        Kk              = PkM*Hk'*pinv(Sk);
        xk              = xkM + Kk*(Accel-hk);
    end    
    Pk              = (eye(9) - Kk*Hk)*PkM;
%     Updating Loops
    eulerEFK(:,ttfwd)     = xk(1:3);
    err(:,ttfwd)       = (Accel-getQT(xk(1:3))*(gN+xk(7:end)));
    PkEKF(:,:,ttfwd)   = Pk; 
    xkEKF(:,ttfwd)     = xk;
    Pkm1            = Pk;
    xkm1            = xk;
end
toc
% Backward
xkS(:,NN)   = xkEKF(:,NN);
% xkS(3,NN)   = imu.realeulerrad(3,NN);
% xkS(7:end,NN)   = gN;

% xkS(:,NN)   = [imu.realeulerrad(:,NN);xkEKF(4:end,NN)];
PkS(:,:,NN)   = PkEKF(:,:,NN);
% keyboard
tic
for ttbkwd = (NN-1):-1:1
    T = dT(ttbkwd);
    PkFwd = PkEKF(:,:,ttbkwd);
    xkFwd = xkEKF(:,ttbkwd);
    [Fk,xkp1M,Bk]     = systemDynamics(xkEKF(:,ttbkwd), Gyro(:,ttbkwd), T, Ca);
    Qksmall              = [Bk*(Qgyro+Qbias)*Bk', -Bk*Qbias;-Qbias*Bk', Qbias];
%     theQT'*Qacc*theQT can be Qdyn
    Qk = [Qksmall, zeros(6,3); zeros(3,6), Qdyn];
    Pkp1M           = Fk*PkFwd*Fk' + Qk;
    Gk              = PkFwd*Fk'*pinv(Pkp1M);
    xkS(:,ttbkwd)     = xkFwd + Gk*(xkS(:,ttbkwd+1)- xkp1M);
    PkS(:,:,ttbkwd)     = PkFwd + Gk*(PkS(:,:,ttbkwd+1)- Pkp1M)*Gk';
    eulerS(:,ttbkwd)    = xkS(1:3,ttbkwd);
    aD(:,ttbkwd)        = xkS(7:end,ttbkwd);
end


toc
    
    figure('name', strcat('Ca=',num2str(Ca),' QaccAmp=',num2str(QaccAmp)));
    subplot(2,1,1);
    plot(rad2deg(imu.realeulerrad'),'--');
    hold on
%     plot(rad2deg(eulerEFK'),'-.','Linewidth',2);
    plot(rad2deg(eulerS'),'Linewidth',2);
    title('Filtered')
%     legend('R True','P True','Y True',...
%         'R Fwd-Bkwd','P Fwd-Bkwd','Y Fwd-Bkwd');
%     results.Bias = omegaBias;
    results.euler = eulerS;
    subplot(2,1,2);
    plot(aD')
    hold on;
    plot(imu.dynaccGlobal','--');
    title('Dynamic Acceleration');
end
%%
%%
function [Fk,xkp1,Bk] = systemDynamics(xk, Omegak, T, Ca)
% xk = Theta = [phi;theta;psi] = [roll;pitch;yaw];
phi = xk(1); theta = xk(2); psi = xk(3);
biask = xk(4:6);
sPh = sin(phi); cPh = cos(phi);
tTh = tan(theta); scTh = 1/cos(theta);
Bk      = [1, sPh*tTh, cPh*tTh;
           0, cPh, -sPh;
           0, sPh*scTh, cPh*scTh];
% parital diffs
Bk_phi  = [0, cPh*tTh, -sPh*tTh;
           0, -sPh, -cPh;
           0, cPh*scTh, -sPh*scTh];
% % % % 
Bk_theta= [0, sPh*scTh^2, cPh*scTh^2;
           0, 0, 0;
           0, sPh*scTh*tTh, cPh*scTh*tTh];
% 
Bk_psi  = zeros(3,3);
unbiasedOmegak = Omegak - biask;
% Fksmall      = eye(6) + [Bk_phi*unbiasedOmegak, Bk_theta*unbiasedOmegak,...
%     Bk_psi*unbiasedOmegak, - Bk; zeros(3,6)]*T;
Fk = [eye(3,3) + [Bk_phi*unbiasedOmegak*T, Bk_theta*unbiasedOmegak*T,...
    Bk_psi*unbiasedOmegak*T], - Bk*T, zeros(3,3);...
    zeros(3,3), eye(3,3), zeros(3,3);...
    zeros(3,6), Ca*eye(3,3)];
xkp1    = [xk(1:6)+[Bk*unbiasedOmegak;zeros(3,1)]*T;Ca*xk(7:end)];
end

%%
function [QT] = getQT(Theta)
    phi = Theta(1); theta = Theta(2); psi = Theta(3);
    Rz_yaw      = [cos(psi), sin(psi), 0;
                   -sin(psi), cos(psi), 0;
                   0 ,0, 1];
    Ry_pitch    = [cos(theta), 0 ,-sin(theta);
                   0, 1,0;
                   sin(theta), 0, cos(theta)];
    Rx_roll     = [1, 0, 0;
                    0, cos(phi), sin(phi);
                    0, -sin(phi), cos(phi)];
    QT              = Rx_roll*Ry_pitch*Rz_yaw;
end

%% 
% gN = gravity in inertial coordinate system => 3X1
% hk = 3x1, Jh = 3x9
function [hk, Jh] = observationDynamicsEKF(xk, gN)
    phi = xk(1); theta = xk(2); psi = xk(3);
    
    Rz_yaw      = [cos(psi), sin(psi), 0;
                   -sin(psi), cos(psi), 0;
                   0 ,0, 1];
    Ry_pitch    = [cos(theta), 0 ,-sin(theta);
                   0, 1,0;
                   sin(theta), 0, cos(theta)];
    Rx_roll     = [1, 0, 0;
                    0, cos(phi), sin(phi);
                    0, -sin(phi), cos(phi)];
%     Rates/Derivatives
    Rz_yaw_rate     = [-sin(psi), cos(psi), 0;
                       -cos(psi), -sin(psi), 0;
                       0 ,0, 0];
    Ry_pitch_rate   = [-sin(theta), 0 ,-cos(theta);
                       0, 0,0;
                       cos(theta), 0, -sin(theta)];
    Rx_roll_rate    = [0, 0, 0;
                       0, -sin(phi), cos(phi);
                       0, -cos(phi), -sin(phi)];
    QT              = Rx_roll*Ry_pitch*Rz_yaw;
    QT_roll         = Rx_roll_rate*Ry_pitch*Rz_yaw;
    QT_pitch        = Rx_roll*Ry_pitch_rate*Rz_yaw;
    QT_yaw          = Rx_roll*Ry_pitch*Rz_yaw_rate;
%     keyboard
    acc_Global      = gN + xk(7:end);
%     Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT];
    Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT];
    hk              = QT*acc_Global;
end
%% 
% gN = gravity in inertial coordinate system => 3X1
% hk = 3x1, Jh = 3x9
function [hk, Jh] = observationDynamicsEKF3(xk, gN)
    phi = xk(1); theta = xk(2); psi = xk(3);
    
    Rz_yaw      = [cos(psi), sin(psi), 0;
                   -sin(psi), cos(psi), 0;
                   0 ,0, 1];
    Ry_pitch    = [cos(theta), 0 ,-sin(theta);
                   0, 1,0;
                   sin(theta), 0, cos(theta)];
    Rx_roll     = [1, 0, 0;
                    0, cos(phi), sin(phi);
                    0, -sin(phi), cos(phi)];
%     Rates/Derivatives
    Rz_yaw_rate     = [-sin(psi), cos(psi), 0;
                       -cos(psi), -sin(psi), 0;
                       0 ,0, 0];
    Ry_pitch_rate   = [-sin(theta), 0 ,-cos(theta);
                       0, 0,0;
                       cos(theta), 0, -sin(theta)];
    Rx_roll_rate    = [0, 0, 0;
                       0, -sin(phi), cos(phi);
                       0, -cos(phi), -sin(phi)];
    QT              = Rx_roll*Ry_pitch*Rz_yaw;
    QT_roll         = Rx_roll_rate*Ry_pitch*Rz_yaw;
    QT_pitch        = Rx_roll*Ry_pitch_rate*Rz_yaw;
    QT_yaw          = Rx_roll*Ry_pitch*Rz_yaw_rate;
%     keyboard
    acc_Global      = gN + xk(7:end);
%     Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT];
%     Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT];
%     hk              = QT*acc_Global;
    Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT;...
                        0,1,0,zeros(1,6)];
    hk              = [QT*acc_Global;theta];
end
%% 
% gN = gravity in inertial coordinate system => 3X1
% hk = 3x1, Jh = 3x9
function [hk, Jh] = observationDynamicsEKF2(xk, gN)
    phi = xk(1); theta = xk(2); psi = xk(3);
    
    Rz_yaw      = [cos(psi), sin(psi), 0;
                   -sin(psi), cos(psi), 0;
                   0 ,0, 1];
    Ry_pitch    = [cos(theta), 0 ,-sin(theta);
                   0, 1,0;
                   sin(theta), 0, cos(theta)];
    Rx_roll     = [1, 0, 0;
                    0, cos(phi), sin(phi);
                    0, -sin(phi), cos(phi)];
%     Rates/Derivatives
    Rz_yaw_rate     = [-sin(psi), cos(psi), 0;
                       -cos(psi), -sin(psi), 0;
                       0 ,0, 0];
    Ry_pitch_rate   = [-sin(theta), 0 ,-cos(theta);
                       0, 0,0;
                       cos(theta), 0, -sin(theta)];
    Rx_roll_rate    = [0, 0, 0;
                       0, -sin(phi), cos(phi);
                       0, -cos(phi), -sin(phi)];
    QT              = Rx_roll*Ry_pitch*Rz_yaw;
    QT_roll         = Rx_roll_rate*Ry_pitch*Rz_yaw;
    QT_pitch        = Rx_roll*Ry_pitch_rate*Rz_yaw;
    QT_yaw          = Rx_roll*Ry_pitch*Rz_yaw_rate;
    acc_Global      = gN + xk(7:end);
    Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT;...
                        0,0,1,zeros(1,6)];
    hk              = [QT*acc_Global;psi];
end