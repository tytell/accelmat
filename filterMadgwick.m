function results = filterMadgwick(imu, varargin)
    %     Options - default
    opt.type            = 'complementary';
    for i=1:2:(nargin-1)
        if strcmp(varargin{i},'alphabybeta'), opt.alphabybeta=varargin{i+1};
        elseif strcmp(varargin{i},'twindow'), opt.twindow=varargin{i+1};
        elseif strcmp(varargin{i},'zeta'), opt.zeta=varargin{i+1};
        elseif strcmp(varargin{i},'type'), opt.type=varargin{i+1}; %options - kalman, complementary
        else error('Invalid argument');
        end
    end
%     Tuning parameters for respective filter type
    if(strcmp(opt.type,'complementary'))
        opt.typenum         = 1;
        opt.alphabybeta     = 10;
        opt.twindow         = 125;
        opt.zeta            = 0.95;
    elseif (strcmp(opt.type,'kalman'))
        opt.typenum     = 2;
        Pkm1            = randn(6,6);
        xkm1Est         = 0.5*rand(6,1);
        samplefreq      = 1e-3;
        gyronoisestd    = deg2rad(0.01*sqrt(0.5/samplefreq)); % unit = radians
        accelnoisestd   = (3e-4*sqrt(0.5/samplefreq));
        Q               = diag([ones(3,1);0.5*ones(3,1)])*gyronoisestd;
        R               = accelnoisestd*eye(3);
    elseif (strcmp(opt.type,'EKF'))
        opt.typenum     = 3;
        Pkm1            = randn(6,6);
        xkm1Est         = 0.5*rand(6,1);
        samplefreq      = 1e-3;
        T               = samplefreq;
        gyronoisestd    = deg2rad(0.01*sqrt(0.5/samplefreq)); % unit = radians
        accelnoisestd   = (3e-4*sqrt(0.5/samplefreq));
        Q               = diag([0.75*ones(3,1);0.5*ones(3,1)])*gyronoisestd;
        R               = accelnoisestd*eye(3);
    else error('Invalid filtering type - input kalman or complemetary');
    end
    
    time            = imu.t;
    Gyroscope       = imu.gyro;
    Accelerometer   = imu.acc;
    DynamicAcceleration_Sensor = zeros(length(time), 3);

    omegaBias = zeros(length(time)+1,3);
    qEstimate = zeros(length(time)+1,4);
    qEstimate(1,:) = [1 0 0 0];
    dT          = [diff(time)]; dT = [dT(1);dT];
    
    tic
    for t = 1:length(time)
%         Step 1 - Acceleromter processing
        Gravity_Sensor   = Accelerometer(t,:) - DynamicAcceleration_Sensor(t,:);
        qAccEst = getQuartGravity_Optimization(Gravity_Sensor, qEstimate(t,:),'constrained');
%         keyboard
%         Step 2 - Gyroscope processing
        if(opt.typenum==1)
%             OPTION 1 - COMPLEMENTARY FILTER
            if(t>opt.twindow)
                gyro_unbiased   = Gyroscope(t,:) - ...
                    opt.zeta*mean(omegaBias((t-opt.twindow):t,:),1);
            else
                gyro_unbiased  = Gyroscope(t,:);
            end
            qDotGyro        = 0.5*quaternProd(qEstimate(t,:),...
                [0, gyro_unbiased]);
            qGyroEst        = qEstimate(t,:) + qDotGyro*dT(t);
    %         Step 3 - find optimal estiamte
    %           fusion=> qEst = gamma*qAccEst + (1-gamma)*qGyroEst
    % %         Variable weight - high omega => use 
            gamma = 1/(1 + opt.alphabybeta*norm(gyro_unbiased));
            qEst = gamma*qAccEst + (1-gamma)*qGyroEst;
    %         Unconstrainted optimization and normalization of the quaternion
            
% % % %              options = optimoptions('fminunc','GradObj','on','Display','notify');
% % % %             qEst = fminunc(@(q)minErrfunc(q,qAccEst, qGyroEst, (1-gamma)/gamma),qEstimate(t,:),options);
            qEst = qEst/norm(qEst);
%         Step 4 - Updating all the stuff
        qEstimate(t+1,:) = qEst;
        euler(t,:) = quatern2euler(quaternConj(qEst));
        qDotEst         = (qEstimate(t+1,:) - qEstimate(t,:))/dT(t);
        qOmegaEst       = 2*quaternProd(quaternConj(qEst),...
            qDotEst);
        omegaBias(t+1,:) = Gyroscope(t,:) - qOmegaEst(2:end);
        qGravity_Sensor = quaternProd(quaternConj(qEst),...
            quaternProd([0 0 0 1],qEst));
        DynamicAcceleration_Sensor(t+1,:) = Accelerometer(t,:) - ...
            qGravity_Sensor(2:end);
        
        elseif(opt.typenum == 2)
%             OPTION 2 - KALMAN FILTER
%             Model => x= [phi;theta;psi; bias (3x3)]
%             x(k+1)  = A x(k) + B u(k) + w(k)
%             z(k)    = H x(k) + v(k)
            phi = xkm1Est(1); theta = xkm1Est(2); psi = xkm1Est(3);
            Hk  = [eye(3,3), zeros(3,3)];
            uk  = Gyroscope(t,:)';
            zk  = quatern2euler(quaternConj(qAccEst))';
%             Kalman Step 1 - Marginal gaussian
            [xkmhat, Jfk] = systemDynamicsEKF(xkm1Est, uk, dT(t));
            Pkm         = Jfk*Pkm1*Jfk' + Q;
%             Kalman Step 2 - Conditional gaussian
            Kk          = Pkm*Hk'*pinv(Hk*Pkm*Hk' + R);
            xkEst       = xkmhat + Kk*(zk-Hk*xkmhat);
            Pk          = (1-Kk*Hk)*Pkm;
            euler(t,:)       = xkEst(1:3);
%             Update for next loop and dynamic acceleration
            xkm1Est = xkEst;
            Pkm1    = Pk;
            R = euler2rotMat(xkEst(1), xkEst(2), xkEst(3));
            DynamicAcceleration_Sensor(t+1,:) = Accelerometer(t,:)- [0,0,1]*R';
        elseif(opt.typenum == 3)
%             OPTION 3 - EXTENDED KALMAN FILTER
%             Model => x= [phi;theta;psi; bias (3x3)]
%             x(k+1)  = F(x(k),u(k)) + w(k)
%             z(k)    = H(x(k)) + v(k)

%             Extended Kalman Step 1 - Marginal gaussian
            Gyro_vector = Gyroscope(t,:)';
            zk = Gravity_Sensor';
            xk = xkm1Est;
            [xkmhat, Jfk] = systemDynamicsEKF(xk, Gyro_vector, T);
            Pkm         = Jfk*Pkm1*Jfk' + Q;
%             Extended Kalman Step 2 - Conditional gaussian
            gN = [0,0,1]';
            [hk, Jh] = observationDynamicsEKF(xkmhat, gN);
            Kk          = Pkm*Jh'*pinv(Jh*Pkm*Jh' + R);
            Pk          = (1-Kk*Jh)*Pkm;
            xkEst       = xkmhat + Kk*(zk-hk);
            euler(t,:)       = xkEst(1:3);
%             Update for next loop and dynamic acceleration
            xkm1Est = xkEst;
            Pkm1    = Pk;
            QR = sensor2inertialRotation(xkEst(1), xkEst(2), xkEst(3));
            DynamicAcceleration_Sensor(t+1,:) = Accelerometer(t,:)- gN'*QR;
        end
    end
    toc
    figure()
    plot(qEstimate);
    plot(rad2deg(imu.realeulerrad'),'--');
    hold on
    plot(rad2deg(euler));
    results.Bias = omegaBias;
    results.euler = euler;
    results.qEstimate = qEstimate;
end

% Outputs the skew symmetric cross product matrix associated with cross
% product operation for a x b => [ax] b where [ax] is a skew symmetric
% matrix. Refer - http://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
% 
% [theMatrix] = crossProductMatrix(theVector)
% Input = theVector - N X 3
% Output = theMatrix - 3 X 3 X N
function [theMatrix] = crossProductMatrix(theVector)
    theMatrix = zeros(3, 3,N);
    for ii = 1:size(theVector,1)
        theMartix(:,:,N) = [0, -theVector(ii,3), theVector(ii,2);
            theVector(ii,3), 0, -theVector(ii,1);
            -theVector(ii,2), theVector(ii,1), 0];
    end
end


%%
function qAccEst = getQuartGravity_Optimization(Gravity_Sensor, q0, varargin)
    if(strcmp(varargin{1},'constrained'))
%         Constrained optimization - ||q|| = 1
        options = optimoptions('fmincon','GradObj','on','display','notify','Algorithm','active-set');
        qConst = fmincon(@(q)costfunc(q,Gravity_Sensor),q0,[],[],[],[],[],[],@theconstraints,options);    
        qAccEst = qConst/norm(qConst);
    else
%         Unconstrainted optimization and normalization of the quaternion
        options = optimoptions('fminunc','GradObj','on','Display','notify');
        qConst = fminunc(@(q)costfunc(q,Gravity_Sensor),q0,options);
        qAccEst = qConst/norm(qConst);
    end
end

%% 
function [f,g] = costfunc(q, Accelerometer)
    F = [2*(q(2)*q(4) - q(1)*q(3)) - Accelerometer(1)
        2*(q(1)*q(2) + q(3)*q(4)) - Accelerometer(2)
        2*(0.5 - q(2)^2 - q(3)^2) - Accelerometer(3)];
    J = [-2*q(3),	2*q(4),    -2*q(1),	2*q(2)
        2*q(2),     2*q(1),     2*q(4),	2*q(3)
        0,         -4*q(2),    -4*q(3),	0];
    f = 0.5*(F'*F);
    g = J'*F;
end
function [c, ceq, GC, GCeq] = theconstraints(q)
    c = [];
    GC = [];
    ceq  = q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2 - 1;
    GCeq = 2*diag(q);
end

% 
function [f,g] = minErrfunc(q,qA,qG,lambda)
    F   = [q-qA,lambda*(q-qG)]';
    J   = [eye(4);lambda*eye(4)]; 
    f = 0.5*(F'*F);
    g   = J'*F;    
end

% This function outputs the mapping matrix between the Roll-Pitch-Yaw rate
% and sensor observed angular rates (i.e. ideal gyroscope observations)
% Theta_dot = B*Omega_sensor
% Theta_dot = [phi_dot, theta_dot, psi_dot]' = [roll_rate, pitch_rate,
% yaw_rate]'
% B = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);]
% Remember : Rotation Matrix = R()
% 
function [B] = sensor2RPYrate(roll, pitch, yaw)
phi = roll; theta = pitch; psi = yaw;
B = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
    0, cos(phi), -sin(phi);
    0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
end

% {e} = Rx(roll)*Ry(pitch)*Rz(yaw)*{E} = Q'*{E}
% {e},{E} => 3x1, Rx,Ry,Rz => 3x3
% However, all vectors are written in vector form e.g. a_sensor = 3x1
% a = {e}'*a_sensor = {E}'*a_inertial
% => a_inertial = Q*a_sensor and Q = Rz(yaw)'*Ry(pitch)'*Rx(roll)'
% =========================
% v_inertial(3x1) = Q*v_sensor(3x1)
% =========================
function [Q] = sensor2inertialRotation(roll,pitch,yaw)
    phi = roll; theta = pitch; psi = yaw;
    Rz_yaw      = [cos(psi), sin(psi), 0;
                   -sin(psi), cos(psi), 0;
                   0 ,0, 1];
    Ry_pitch    = [cos(theta), 0 ,-sin(theta);
                   0, 1,0;
                   sin(theta), 0, cos(theta)];
    Rx_roll     = [1, 0, 0;
                    0, cos(phi), sin(phi);
                    0, -sin(phi), cos(phi)];
    Q           = Rz_yaw'*Ry_pitch'*Rx_roll';
end

% xk = 6x1
% Jf = 6x6
function [xkp1,Jf] = systemDynamicsEKF(xk, Gyro, T)
%     x = [Theta;bias]
    Gyro_unbiased   = Gyro - xk(4:end); 
    phi = xk(1); theta = xk(2); psi = xk(3);
    B = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
    0, cos(phi), -sin(phi);
    0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
% % % % % 
    B_x1 = [0, cos(phi)*tan(theta), -sin(phi)*tan(theta);
    0, -sin(phi), -cos(phi);
    0, cos(phi)*sec(theta), -sin(phi)*sec(theta)];
% % % % 
    B_x2 = [0, sin(phi)*(sec(theta)^2), cos(phi)*(sec(theta)^2);
    0, 0, 0;
    0, sin(phi)*sec(theta)*tan(theta), cos(phi)*sec(theta)*tan(theta)];
% 
    Jf = eye(6) + [B_x1*Gyro_unbiased, B_x2*Gyro_unbiased, zeros(3,1), -B;
        zeros(3,6)]*T;
    xkp1 = xk + [B*Gyro_unbiased; zeros(3,1)];
end

% gN = gravity in inertial coordinate system => 3X1
% hk = 3x1, Jh = 3x6
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
    Rz_yaw_rate     = [cos(psi), sin(psi), 0;
                       -sin(psi), cos(psi), 0;
                       0 ,0, 1];
    Ry_pitch_rate   = [cos(theta), 0 ,-sin(theta);
                       0, 1,0;
                       sin(theta), 0, cos(theta)];
    Rx_roll_rate    = [1, 0, 0;
                       0, cos(phi), sin(phi);
                       0, -sin(phi), cos(phi)];
    QT              = Rx_roll*Ry_pitch*Rz_yaw;
    QT_roll         = Rx_roll_rate*Ry_pitch*Rz_yaw;
    QT_pitch        = Rx_roll*Ry_pitch_rate*Rz_yaw;
    QT_yaw          = Rx_roll*Ry_pitch*Rz_yaw_rate;
    Jh              = [QT_roll*gN, QT_pitch*gN, QT_yaw*gN, zeros(3,3)];
    hk              = QT*gN;
end