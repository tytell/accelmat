function imu = get_orient_imu(imu, varargin)

opt.method = 'simple';
opt.gyrooffset = [-16 -8];
opt.gyroband = [0.1 10];
opt.checkacc = false;

opt = parsevarargin(opt,varargin,2);

imu.rate = 1/mean(diff(imu.t));

[b,a] = butter(5,opt.gyroband(1)/(imu.rate/2),'low');
gyrolo = filtfilt(b,a, imu.gyro);

gyros = imu.gyro - gyrolo;
[b,a] = butter(5,opt.gyroband(2)/(imu.rate/2), 'low');
gyros = filtfilt(b,a, gyros);

orient = cumtrapz(imu.t, gyros);
orient(:,1) = orient(:,1) - opt.gyrooffset(:,1);
orient(:,2) = orient(:,2) - opt.gyrooffset(:,2);

if opt.checkacc
    [b,a] = butter(5,0.5/(imu.rate/2),'low');
    acclo = filtfilt(b,a,imu.acc);
    mag = sqrt(sum(acclo.^2,2));
    acclo = acclo ./ repmat(mag,[1 3]);

    orienta(:,1) = (pi + unwrap(atan2(acclo(:,2),acclo(:,3)))) * 180/pi;
    if (any(orienta(:,1) > 360))
        orienta(:,1) = orienta(:,1) - 360;
    end
    orienta(:,2) = -(pi + unwrap(atan2(acclo(:,1),acclo(:,3)))) * 180/pi;
    
    plot(imu.t, orient(:,1:2));
    addplot(imu.t, orienta);
end

orient = orient * pi/180;

gvec = [sin(orient(:,2)) -sin(orient(:,1)) -cos(orient(:,2))];
gvecmag = sqrt(sum(gvec.^2,2));
gvec = gvec ./ repmat(gvecmag,[1 3]);

imu.orient = orient;
imu.acchi = imu.acc - gvec;
