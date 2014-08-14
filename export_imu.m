function export_imu(filename,calibfile, outfile, varargin)

opt.resamplerate = 250;
opt.kinematicsfile = '';
opt = parsevarargin(opt,varargin, 4);

if nargin == 0
    [fn,pn] = uigetfile('*.h5','Choose IMU data file');
    filename = fullfile(pn,fn);
    
    [fn,pn] = uigetfile('*.mat','Choose calibration file',pn);
    calibfile = fullfile(pn,fn);
    
    [fn,pn] = uigetfile('*.mat','Choose kinematics file (or Cancel for none)',pn);
    if ~isempty(fn)
        kinfile = fullfile(pn,fn);
    else
        kinfile = '';
    end
    
    [fn,pn] = uiputfile('*.txt','Choose output file',pn);
    outfile = fullfile(pn,fn);
else
    kinfile = opt.kinematicsfile;
end;

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
X0(:,4:6) = X0(:,4:6) * 9.800;        %convert to m/s^2
for i = 1:size(X0,2)
    X1(:,i) = interp1(t0,X0(:,i), t1);
end
X1 = [t1 X1];

chantitle = {'Pitch','Roll','Yaw','UpAcc','FwdAcc','SideAcc'};
tplt = '%.5f\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.5f';
if ~isempty(kinfile)
    K = load(kinfile,'haxmmss','haymmss','humms','hvmms','hxmm',...
        'hymm','txmm','tymm','t');
    
    X2 = [K.hxmm' K.hymm' K.txmm' K.tymm' K.humms' K.hvmms' K.haxmmss' K.haymmss'];
    K.t = K.t(end) - K.t;
    for i = 1:size(X2,2)
        X3(:,i) = interp1(K.t,X2(:,i), t1);
    end
    X3 = X3 / 1000;       % convert to m
    
    X1 = [X1 X3];
    
    chantitle = [chantitle, {'FwdPos','SidePos','TailFwd','TailSide','FwdVel','SideVel',...
        'FwdAccVid','SideAccVid'}];
    tplt = [tplt '\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f'];
end

fid = fopen(outfile,'w');
titletplt = ['ChannelTitle=\t' repmat('%s\t',[1 length(chantitle)-1]) '%s\n'];
fprintf(fid, titletplt, chantitle{:});
fprintf(fid, [tplt '\n'], X1');
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
if exist('X3','var')
    fwdpos = X3(:,1);
    sidepos = X3(:,2);
    tailfwd = X3(:,3);
    tailside = X3(:,4);
    fwdvel = X3(:,5);
    sidevel = X3(:,6);
    fwdaccvid = X3(:,7);
    sideaccvid = X3(:,8);
    
    extravars = {'fwdpos','sidepos','tailfwd','tailside','fwdvel','sidevel',...
        'fwdaccvid','sideaccvid'};
else
    extravars = {};
end
save(outmatfile,'t','pitch','roll','yaw','accup','accfwd','accside', extravars{:});




    