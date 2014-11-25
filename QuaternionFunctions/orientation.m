function [qrotV qrotA]  = orientation(imu, Gimu)
% imu: use load_imu or load_2imu to get acceleration and gyro data
% grav: measured acceleration of the chip at rest (measure of gravity)

% Initiate variables
    grav = mean(Gimu.acc);
    errA = max(std(Gimu.acc));
    angAcc = mean(Gimu.gyro);
    errOmega = max(std(Gimu.gyro));
    a0 = grav/norm(grav);
%   Z = a0/(norm(a0));
    Z = [0 0 1];
    
    A0norm = repmat(sqrt(sum(Gimu.acc.*Gimu.acc,2)),1,3);
    A0unit = Gimu.acc./A0norm;
    Theta0 = acos(A0unit(:,3));
    errTheta = std(Theta0);
    
    % was getting a number higher then 1 for the dot product below.
    % assumed measurement was in degrees and needed to be converted to
    % radians
    % theta0 = acos((dot(a0, Z)*pi)/180);
    theta0 = acos(dot(a0, Z));
    V0 = cross(-a0, Z);
    P = sin(theta0/2);
    Q = [V0/mag(V0)];
    q(1,:) = [cos(theta0/2), P*Q];
    f = 0.01;
    qrotA = [0 0 0];
 % Generate Quaternions
    for i = 2:length(imu.t)
        ig = imu.gyro(i,1); ia = imu.acc(i,1);
        jg = imu.gyro(i,2); ja = imu.acc(i,2);
        kg = imu.gyro(i,3); ka = imu.acc(i,3);
%        P = imu.gyro(i-1,:);
%        Q = imu.gyro(i,:);
        AngV = [ig,jg,kg]; Acc = [ig,jg,kg];
        theta = acos((dot(Acc, Z)*pi)/180);
        Vect = Acc;
        % want us to multiply a quaternion by the angular velocity but that
        % is only a 1x3 vector and quatmultiply requires two 1x4
        % quaternions. Not sure what to do here
        QQ = quatmultiply2(q(i-1,:),[0,(AngV/f)]);
        Omega = quatmultiply2(QQ,quatinv2(q(i-1,:)));
        Phi = [1, ([Omega(2), Omega(3), Omega(4)]/2)];
        q(i,:) = quatmultiply2(Phi, q(i-1,:));
      
        % METHOD 2 Start --------------------------------------------------
%         if (mag(a0-Acc) < mag(Acc) && mag(Acc) < mag(a0+Acc) && ...
%             acos(abs(kg)) < theta) || ...
%            (mag(a0-Acc) < mag(Acc) && mag(Acc) < mag(a0+Acc) && ...
%             mag(AngV) < AngV)
        if (abs(norm(grav)- norm(Acc)) < errA && acos(abs(Acc(1,3)))-theta0 < errTheta) || ...
           (abs(norm(grav)- norm(Acc)) < errA && abs(norm(angAcc)-norm(AngV)) < errOmega) 
            oldQ = q(i,:);
            % Same here, Acc is only 1x3
            A = quatmultiply2(quatmultiply2(q(i,:), [0, Acc]), quatinv2(q(i,:)));
            Psi = [1, A(3)/2, -A(1)/2, 0];
            q(i,:) = quatmultiply2(Psi, oldQ);
        end
        % METHOD 2 End ----------------------------------------------------

% %        theta = acos(((dot(P,Q))/(mag(P)*mag(Q))));        
% %        u = cos(theta/2) + sin(theta/2) * ((ii+jj+kk)/mag(Q));
% %        [EMx EMy EMz] = EulerMat([ii jj kk]);
        
        qrotV(i,:) = quatrotate2(q(i,:),Vect);
    end

        qrotA = q;
end

