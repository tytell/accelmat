%first set up calibration
[imu,calib] = load_imu('rawdata/Bg5-Calib-001.h5',[],'clickzero','calib');

%then load the imu data
imu = load_imu('rawdata/Bg5-20_0hz-005.h5',calib);

%pull out the appropriate variables
t = imu.t;
acc = imu.acc;
gyro = imu.gyro;

%and save the imu file
save Bg5-20_0hz-005.mat

%then we can sync
sync_video_and_data
