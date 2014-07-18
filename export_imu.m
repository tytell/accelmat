function export_imu(filename,calibfile, outfile, varargin)

opt.resamplerate = 250;
opt = parsevarargin(opt,varargin, 4);

if nargin == 0
    [fn,pn] = uigetfile('*.h5','Choose IMU data file');
    filename = fullfile(pn,fn);
    
    [fn,pn] = uigetfile('*.mat','Choose calibration file',pn);
    calibfile = fullfile(pn,fn);
    
    [fn,pn] = uiputfile('*.txt','Choose output file',pn);
    outfile = fullfile(pn,fn);
end

if (ischar(calibfile) && exist(calibfile,'file'))
    load(calibfile,'calib');
elseif isstruct(calibfile)
    calib = calibfile;
end
imu = load_imu(filename,calib);

imu = get_orient_imu(imu, 'getoffset');

t0 = imu.t;
t1 = (t0(1):1/opt.resamplerate:t0(end))';
X0 = [imu.orient imu.acchi];
X0(:,1:3) = X0(:,1:3) * 180/pi;     %convert to degrees
X0(:,4:6) = X0(:,4:6) * 9.8;        %convert to m/s^2
for i = 1:size(X0,2)
    X1(:,i) = interp1(t0,X0(:,i), t1);
end
X1 = [t1 X1];

fid = fopen(outfile,'w');
fprintf(fid, 'ChannelTitle=\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'Pitch','Roll','Yaw','Up','Forward','Side');
fprintf(fid, '%.5f\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.5f\n', X1');
fclose(fid);

[pn,fn,ext] = fileparts(outfile);
outmatfile = fullfile(pn,[fn '.mat']);

t = t0;
pitch = X0(:,1);
roll = X0(:,2);
yaw = X0(:,3);
accup = X0(:,4);
accfwd = X0(:,5);
accside = X0(:,6);
save(outmatfile,'t','pitch','roll','yaw','accup','accfwd','accside');




    