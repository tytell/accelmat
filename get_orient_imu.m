function imu = get_orient_imu(imu, varargin)

opt.method = 'madgwick';        % or 'simple' or 'ertss'
opt.gyrooffset = [-16 -8];
opt.gyroband = [0.1 10];
opt.accband = [0 30];
opt.getoffset = false;
opt.imuposition = [6.6 11.4 -7];        % distance from COM in mm
opt.quiet = false;

opt.beta = 2.86;            % deg/sec
opt.initwindow = 0.5;       % use the first 0.5 sec to initialize orientation

opt.Qgyro = [];
opt.constbiasgyro = [];
opt.Qacc = [];
opt.Qdyn = [];
opt.Qbias = [];
opt.Ca = 1.1;
opt.knownorientind = [];

opt = parsevarargin(opt,varargin,2);

imu.rate = 1/mean(diff(imu.t));
dt = 1/imu.rate;

N = length(imu.t);

%ideally we want a bandpass, but matlab doesn't do well with very
%low cutoffs, so we do a running mean type operation to get rid of
%the very low frequencies
if (opt.gyroband(1) > 0)
    gyrolo = get_low_baseline(imu.t, imu.gyro, opt.gyroband(1));
    
    gyros = imu.gyro - gyrolo;
else
    gyros = imu.gyro;
end
%and then a low pass filter to get rid of the high frequencies
[b,a] = butter(5,opt.gyroband(2)/(imu.rate/2), 'low');
gyros = filtfilt(b,a, gyros);
gyros = gyros - repmat(nanmean(gyros),[size(imu.gyro,1) 1]);

switch lower(opt.method)        
    case {'simple','madgwick'}
        if strcmp(opt.method,'simple')
            opt.beta = 0;
        end
        
        %have to be in radians
        gyros = gyros*pi/180;
        opt.beta = opt.beta * pi/180;
        
        %filter the accelerometer data
        if (opt.accband(1) > 0)
            acclo = get_low_baseline(imu.t, imu.acc, opt.accband(1));
            
            accs = imu.acc - acclo;
        else
            accs = imu.acc;
        end
        
        [b,a] = butter(5,opt.accband(2)/(imu.rate/2), 'low');
        accs = filtfilt(b,a, accs);
        
        qorient = zeros(N,4);
        
        %get the initial quaternion, assuming that the first half second of
        %accelerometer data is a correct estimate of gravity
        %flip the sign, because we want our orientation quaternion to tell
        %us the direction of up
        isfirst = imu.t <= imu.t(1) + opt.initwindow;
        acc1 = nanmean(accs(isfirst,:));
        acc1 = -acc1 / norm(acc1);
        ax = acc1(1);
        ay = acc1(2);
        az = acc1(3);
        
        AXZ = ax*sqrt(1 - az);
        AXY = sqrt(ax^2 + ay^2);
        q0 = [0 ...
            AXZ/(sqrt(2)*AXY) ...
            ay*AXZ/(sqrt(2)*ax*AXY) ...
            ax*AXY / (sqrt(2)*AXZ)];
        
        qorient(1,:) = q0 / norm(q0);
        
        qgyro = zeros(N,4);
        gvec = zeros(N,3);
        
        for i = 2:N
            qprev = qorient(i-1,:);
            
            %flip the sign here
            acc1 = accs(i,:);
            acc1 = -acc1 / norm(acc1);

            %quaternion angular change from the gyro
            qdotgyro = 0.5 * quaternProd(qprev, [0 gyros(i,:)]);
            qgyro(i,:) = qprev + qdotgyro * dt;
            
            % Gradient descent algorithm corrective step
            F = [2*(qprev(2)*qprev(4) - qprev(1)*qprev(3)) - acc1(1)
                2*(qprev(1)*qprev(2) + qprev(3)*qprev(4)) - acc1(2)
                2*(0.5 - qprev(2)^2 - qprev(3)^2) - acc1(3)];
            J = [-2*qprev(3),	2*qprev(4),    -2*qprev(1),	2*qprev(2)
                2*qprev(2),     2*qprev(1),     2*qprev(4),	2*qprev(3)
                0,         -4*qprev(2),    -4*qprev(3),	0    ];
            step = (J'*F);
            step = step / norm(step);	% normalise step magnitude
            
            qdot = qdotgyro - opt.beta * step';
            
            qorient(i,:) = qprev + qdot * dt;
            qorient(i,:) = qorient(i,:) / norm(qorient(i,:));
            
            %get the gravity vector from the orientation
            %remember gravity is -Z
            g1 = quaternProd(quaternConj(qorient(i,:)), quaternProd([0 0 0 -1], qorient(i,:)));
            gvec(i,:) = g1(2:4);
        end
        
        imu.orient = quatern2euler(quaternConj(qorient)) * 180/pi;
        imu.rotmat = quatern2rotMat(quaternConj(qorient));
        imu.qorient = qorient;
        imu.gvec = gvec;
        imu.accdyn = accs - gvec;
        
    case 'ertss'
        %have to be in radians
        gyros = gyros*pi/180;
        
        %filter the accelerometer data
        if (opt.accband(1) > 0)
            acclo = get_low_baseline(imu.t, imu.acc, opt.accband(1));
            
            accs = imu.acc - acclo;
        else
            accs = imu.acc;
        end
        
        [b,a] = butter(5,opt.accband(2)/(imu.rate/2), 'low');
        accs = filtfilt(b,a, accs);

        if isempty(opt.Qgyro) || isempty(opt.Qacc)
            error('Must have values for gyro and acceleration covariances (''Qgryo'' and ''Qacc'')');
        end
        
        if isempty(opt.Qdyn)
            opt.Qdyn            = 10*opt.Qacc;
        end
        if isempty(opt.Qbias)
            opt.Qbias           = 0.05*opt.Qgyro;
        end
        
        imu1 = imu;
        imu1.gyro = gyros;
        imu1.acc = accs;
        
        imu = get_orient_imu_ERTSS(imu1, opt.Qgyro,opt.Qbias,opt.Qacc,opt.Qdyn,opt.Ca, 'quiet',opt.quiet, ...
            'knownorientind',opt.knownorientind);
end

angacclo = deriv(imu.t, gyros);

% Convert the units first - deg to radians
imu.gyrocross       = crossProductMatrix(gyros*pi/180);
imu.angacclocross   = crossProductMatrix(angacclo*pi/180);
imu.correctiondist  = opt.imuposition*1e-3;  % Correction distance (1,3) in meters
for ii = 1:size(imu.acc,1)
    imu.accdyn(ii,:)  = imu.accdyn(ii,:) + imu.correctiondist*...
        (imu.angacclocross(:,:,ii) + imu.gyrocross(:,:,ii)^2)';
end

% Outputs the skew symmetric cross product matrix associated with cross
% product operation for a x b => [ax] b where [ax] is a skew symmetric
% matrix. Refer - http://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
% 
% [theMatrix] = crossProductMatrix(theVector)
% Input = theVector - N X 3
% Output = theMatrix - 3 X 3 X N
function [theMatrix] = crossProductMatrix(theVector)

N = size(theVector,1);
theMatrix = zeros(3, 3,N);
for ii = 1:size(theVector,1)
    theMartix(:,:,N) = [0, -theVector(ii,3), theVector(ii,2);
        theVector(ii,3), 0, -theVector(ii,1);
        -theVector(ii,2), theVector(ii,1), 0];
end


