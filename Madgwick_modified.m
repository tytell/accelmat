% Modified Madgwick IMU

function [results] = Madgwick_modified(imu, alphabybeta, twindow, zeta)
%     close all; clear all; clc;
%     clc
%     imu = simulatedData;
%     zeta = 1;
tic
    time            = imu.t;
    Gyroscope       = imu.gyro;
    Accelerometer   = imu.acc;

    DynamicAcceleration_Sensor = zeros(length(time), 3);

    omegaBias = zeros(length(time)+1,3);
    qEstimate = zeros(length(time)+1,4);
    qEstimate(1,:) = [1 0 0 0];
% %     qEstimate(1,:) = randn(4,1);
% %     qEstimate(1,:) = qEstimate(1,:)/norm(qEstimate(1,:));
    dT          = [diff(time)]; dT = [dT(1);dT];
    for t = 1:length(time)
%         Step 1 - Gyroscope processing
        if(t>twindow)
            gyro_unbiased   = Gyroscope(t,:) - zeta*mean(omegaBias((t-twindow):t,:),1);
        else
            gyro_unbiased  = Gyroscope(t,:);
        end
        qDotGyro        = 0.5*quaternProd(qEstimate(t,:),...
            [0, gyro_unbiased]);
        qGyroEst        = qEstimate(t,:) + qDotGyro*dT(t);
%         Step 2 - Acceleromter processing
% % % %         In case there is need for rotation from S(t-1) to S(t)
% % %         qtm12t  = quaternProd(qGyroEst,quaternConj(qEstimate(t,:)));
% % %         modifiedDyAccel = quaternProd(quaternConj(qtm12t), ...
% % %             quaternProd([0, DynamicAcceleration_Sensor(t,:)],qtm12t));
% % %         Gravity_Sensor   = Accelerometer(t,:) - modifiedDyAccel(2:end);
        Gravity_Sensor   = Accelerometer(t,:) - DynamicAcceleration_Sensor(t,:);
        options = optimoptions('fminunc','GradObj','on','Display','notify');
%         Unconstrainted optimization and normalization of the quaternion
        qAccEst = fminunc(@(q)costfunc(q,Gravity_Sensor),qEstimate(t,:),options);
        qAccEst = qAccEst/norm(qAccEst);
% % %         Constrained optimization - ||q||=1
% %         options = optimoptions('fmincon','GradObj','on','display','notify','Algorithm','active-set');
% %         qAccEst = fmincon(@(q)costfunc(q,Gravity_Sensor),qEstimate(t,:),[],[],[],[],[],[],@theconstraints,options);
%         Step 3 - find optimal estiamte
%           fusion=> qEst = gamma*qAccEst + (1-gamma)*qGyroEst
% %         Variable weight
%         gamma = 1/(1+ alphabybeta*norm(qDotGyro));
% %         Constant weight
        gamma = alphabybeta;
%         Unconstrainted optimization and normalization of the quaternion
        qEst = fminunc(@(q)minErrfunc(q,qAccEst, qGyroEst, (1-gamma)/gamma),qEstimate(t,:),options);
        qEst = qEst/norm(qEst);
% % %         Constrained Optimizatio - ||q||=1
% %         options = optimoptions('fmincon','GradObj','on','display','notify','Algorithm','active-set');
% %         qEst = fmincon(@(q)minErrfunc(q,qAccEst, qGyroEst, (1-gamma)/gamma),qEstimate(t,:),[],[],[],[],[],[],@theconstraints,options);
% % %         
% % % %         Another approach - weighted but without optimization
% % %         qEst = gamma*qAccEst + (1-gamma)*qGyroEst;
% % %         qEst = qEst/norm(qEst);
% % % 
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
    end
    toc
    figure('Name', strcat('\alpha / \beta = ',num2str(alphabybeta),' T_{window} = ', num2str(twindow)))
    plot(rad2deg(imu.realeulerrad'),'--');
    hold on
    plot(rad2deg(euler));
    results.Bias = omegaBias;
    results.euler = euler;
    results.qEstimate = qEstimate;
%     keyboard
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