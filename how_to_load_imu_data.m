[~, calib] = load_imu('path/to/calib.h5', [], 'calib', 'clickzero');

% timerange samples relative to trigger, from 5 sec before to 1 sec
% after in the example below
% 
imu = load_imu('path/to/data.h5', calib, ...
    'timerange', [-5 1], ...        % just load from -5 to 1 sec, rel to trigger
    'resamplerate',100);

h(1) = subplot(2,1,1);
plot(imu.t, imu.gyro);

h(2) = subplot(2,1,2);
plot(imu.t, imu.acc);

linkaxes(h, 'x');
