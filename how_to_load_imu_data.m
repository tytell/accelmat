[~, calib] = load_imu('path/to/calib.h5', [], 'calib', 'clickzero');

% timerange samples relative to trigger, from 5 sec before to 1 sec
% after in the example below
% 
imu = load_imu('path/to/data.h5', calib, ...
    'timerange', [-5 1], ...        % just load from -5 to 1 sec, rel to trigger
    'resamplerate',100);

h(1) = subplot(3,1,1);
plot(imu.t, imu.gyro);

h(2) = subplot(3,1,2);
plot(imu.t, imu.acc);


%%% Loading LabChart daat
% Export the relevant part of your data set as a Matlab file from LabChart

importLabChart('/path/to/exported/LabChart/file.mat', ...
    '/path/to/better/matfile.mat');

emg = load('/path/to/better/matfile.mat');

h(3) = subplot(3,1,3);
t0 = 10;        % find the time for the trigger from the EMG data
plot(emg.t - t0, emg.channelname1);     % where channelname1 is the channel name from LabChart

linkaxes(h, 'x');
