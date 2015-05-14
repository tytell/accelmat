function imu = get_orient_imu_ERTSS(imu, Qgyro, Qbias, Qacc, Qdyn, Ca, varargin)

opt.knownorientind = [];        % indices of frames with known orientation
opt.doplot = false;
opt.quiet = false;
opt.initwindow = 0.5;

opt = parsevarargin(opt, varargin, 7);

knownN = opt.knownorientind;

%check for radians vs degrees
if all(abs(imu.gyro(:)) < 250 * pi/180)
    warning('It appears that the gyro readings are already in rad/sec. Results may be strange');
end

% Step 1 - get all data
%convert gyro readings and covariances to rad/sec
Gyro    = imu.gyro' * pi/180; 
Qgyro = Qgyro * pi/180;
Qbias = Qbias * pi/180;

Acc     = imu.acc';
Rk      = Qacc;
xkm1      = zeros(9,1);

%get the initial orientation, assuming that the accelerometer is correct in
%the first few readings
if ~isempty(opt.initwindow)
    isfirst = imu.t <= imu.t(1) + opt.initwindow;
    accmn = nanmean(imu.acc(isfirst,:));
    accmn = accmn / norm(accmn);
    
    r = atan2(accmn(2),accmn(3));
    p = atan2(accmn(1),accmn(3));
    
    xkm1(1:3) = [r p 0];
end

dT       = diff(imu.t);
dT       = [dT(1); dT];
T       = dT(1);
Pkm1      = [(Qgyro+Qbias)*T^2, -Qbias*T, zeros(3,3);-Qbias*T, Qbias,zeros(3,3);...
    zeros(3,6), Qdyn];
gN      = [0,0,-1]';

gN = gN/norm(gN);

NN          = length(Gyro);
eulerEFK       = zeros(3,NN);
eulerS       = zeros(3,NN);
err         = zeros(3,NN);
PkEKF       = zeros(9,9,NN);

i = 1;
progress(1,2*NN,'Forward filter...', 'quiet',opt.quiet);
for ttfwd=1:NN
    T = dT(ttfwd);
    Omegak          = Gyro(:,ttfwd);
    Accel           = Acc(:,ttfwd);
    [Fk,xkM,Bk]     = systemDynamics(xkm1, Omegak, T,Ca);
    theQT           = getQT(xkM(1:3));
    Qk              = [Bk*(Qgyro+Qbias)*Bk'*T^2, -Bk*Qbias*T, zeros(3,3);...
        -Qbias*Bk'*T, Qbias, zeros(3,3);
        zeros(3,6), Qdyn];
    PkM             = Fk*Pkm1*Fk' + Qk;

    if (ismember(ttfwd, knownN))
        %update with known orientation
        [hk2, Jh2] = observationDynamicsEKF2(xkM, gN);
        Hk          = Jh2;
        Sk              = Hk*PkM*Hk' + [Rk,zeros(3,1);zeros(1,4)];
        Kk              = PkM*Hk'*pinv(Sk);
        xk          = xkM + Kk*([Accel;imu.realeulerrad(3,ttfwd)]...
                                -hk2);
    else 
        %update and estimate orientation
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

    progress(i);
    i = i+1;
end

% Backward
xkS(:,NN)   = xkEKF(:,NN);
PkS(:,:,NN)   = PkEKF(:,:,NN);
aD(:,NN)    = xkS(7:end,NN);

for ttbkwd = (NN-1):-1:1
    T = dT(ttbkwd);
    PkFwd = PkEKF(:,:,ttbkwd);
    xkFwd = xkEKF(:,ttbkwd);
    [Fk,xkp1M,Bk]     = systemDynamics(xkEKF(:,ttbkwd), Gyro(:,ttbkwd), T, Ca);
    Qksmall              = [Bk*(Qgyro+Qbias)*Bk', -Bk*Qbias;-Qbias*Bk', Qbias];
    Qk = [Qksmall, zeros(6,3); zeros(3,6), theQT'*Qacc*theQT];
    Pkp1M           = Fk*PkFwd*Fk' + Qk;
    Gk              = PkFwd*Fk'*pinv(Pkp1M);
    xkS(:,ttbkwd)     = xkFwd + Gk*(xkS(:,ttbkwd+1)- xkp1M);
    PkS(:,:,ttbkwd)     = PkFwd + Gk*(PkS(:,:,ttbkwd+1)- Pkp1M)*Gk';
    eulerS(:,ttbkwd)    = xkS(1:3,ttbkwd);
    aD(:,ttbkwd)        = xkS(7:end,ttbkwd);

    progress(i);
    i = i+1;
end

if opt.doplot
    figure('name', strcat('Ca=',num2str(Ca)));
    subplot(2,1,1);
    plot(imu.t,rad2deg(eulerEFK'),'-.','Linewidth',2);
    hold on
    plot(imu.t,rad2deg(eulerS'),'Linewidth',2);
    hold off
    legend('roll','pitch','yaw');
    title('Filtered')
    subplot(2,1,2);
    plot(imu.t,aD')
    title('Dynamic Acceleration');
end

imu.orient = eulerS' * 180/pi;
imu.accdyn = aD';
imu.gyrodrift = xkS(4:6,:)' * 180/pi;


%------------------------------------------------------------------------
function [Fk,xkp1,Bk] = systemDynamics(xk, Omegak, T, Ca)

phi = xk(1); 
theta = xk(2); 
psi = xk(3);

biask = xk(4:6);

sPh = sin(phi); cPh = cos(phi);
tTh = tan(theta); scTh = 1/cos(theta);
Bk      = [1, sPh*tTh, cPh*tTh;
           0, cPh, -sPh;
           0, sPh*scTh, cPh*scTh];
% partial diffs
Bk_phi  = [0, cPh*tTh, -sPh*tTh;
           0, -sPh, -cPh;
           0, cPh*scTh, -sPh*scTh];

Bk_theta= [0, sPh*scTh^2, cPh*scTh^2;
           0, 0, 0;
           0, sPh*scTh*tTh, cPh*scTh*tTh];

Bk_psi  = zeros(3,3);
unbiasedOmegak = Omegak - biask;

Fk = [eye(3,3) + [Bk_phi*unbiasedOmegak*T, Bk_theta*unbiasedOmegak*T,...
    Bk_psi*unbiasedOmegak*T], - Bk*T, zeros(3,3);...
    zeros(3,3), eye(3,3), zeros(3,3);...
    zeros(3,6), Ca*eye(3,3)];
xkp1    = [xk(1:6)+[Bk*unbiasedOmegak;zeros(3,1)]*T;Ca*xk(7:end)];


%------------------------------------------------------------------------
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


%------------------------------------------------------------------------
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

acc_Global      = gN + xk(7:end);

Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT];
hk              = QT*acc_Global;

    % 
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

    acc_Global      = gN + xk(7:end);

    Jh              = [QT_roll*acc_Global, QT_pitch*acc_Global, QT_yaw*acc_Global, zeros(3,3), QT;...
                        0,1,0,zeros(1,6)];
    hk              = [QT*acc_Global;theta];

%------------------------------------------------------------------------
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