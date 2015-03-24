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
        elseif(opt.typenum == 2)
%             OPTION 2 - KALMAN FILTER
%             Model => x= [phi;theta;psi; bias (3x3)]
%             x(k+1)  = A x(k) + B u(k) + w(k)
%             z(k)    = H x(k) + v(k)
            phi = xkm1Est(1); theta = xkm1Est(2); psi = xkm1Est(3);
            W   = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
                0, cos(phi), -sin(phi);
                0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
            Ak  = [eye(3,3), -W*dT(t);
                zeros(3,3) + eye(3,3)];
            Bk  = [W*dT(t);zeros(3,3)];
            Hk  = [ones(1,3), zeros(1,3)];
            uk  = Gyroscope(t,:);
            zk  = qAccEst;
%             Kalman Step 1 - Marginal gaussian
            Pkm         = Ak*Pkm1*Ak' + Q;
            xkmhat      = Ak*xkm1Est + Bk*uk;
%             Kalman Step 2 - Conditional gaussian
            Kk          = Pkm*Hk'*pinv(Hk*Pkm*Hk' + R);
            xkEst       = xkmhat +Kk*(zk-Hk*xkmhat);
            Pk          = (1-Kk*Hk)*Pkm;
            euler       = xkEst(1:3);
%             Update for next loop
            xkm1Est = xkEst;
            Pkm1    = Pk;
        end
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
%     keyboard
    figure()
    plot(qEstimate);
    figure('Name', strcat('alpha/beta = ',num2str(opt.alphabybeta),' T_{window} = ', num2str(opt.twindow)))
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
