function process_bg36_accel_data

accmethod = 'madgwick';
smoothdur = 0.5;
leftsidevortexsign = 1;

fishlen = 161;       % mm

sampfreq = 100;

imuposition = [0.5 10.7 -16.6];        % x y z coords in mm

steadythresh = 0.03;    

acccalibfile = 'F:\Accelerometer h5\accelerometer bg25/calib001.h5'; %'F:\Accelerometer h5\accelerometer bg38/calib001.h5'; %fill in here
massdistfile = 'C:\Code\accelmat\fishmass.mat'; % or 'C:\Code\accelmat\fishmass.mat';
noisefile = 'F:\Accelerometer h5\accelerometer bg25/bg25 and 36 10 min test 001.h5'; % or 'F:\Accelerometer h5\accelerometer bg38/bg38 10min test 001.h5';

baseaccdir = 'F:\Accelerometer h5\accelerometer bg25';    % or 'F:\Accelerometer h5\accelerometer bg38'
basekindir = 'F:\Digitize Fish\digitize fish bg25';    % or 'F:\Digitize Fish\digitize fish bg38'
basepivdir = 'F:\Analyze PIV\analyzepiv bg25';    % or 'F:\Digitize Fish\digitize fish bg38'

accfiles = {

    'Bg25_027.h5'
    'Bg25_028.h5'
    'Bg25_029.h5'
    'Bg25_030.h5'
    'Bg25_031.h5'
    'Bg25_032.h5'
    'Bg25_033.h5'
    'Bg25_034.h5'
    };
kinfiles = { 
    'bg25 28_7Hz 027.mat'
    'bg25 28_7Hz 028.mat'
    'bg25 28_7Hz 029.mat'  
    'bg25 28_7Hz 030.mat'
    'bg25 28_7Hz 031.mat'
    'bg25 28_7Hz 032.mat'
    'bg25 28_7Hz 033.mat'
    'bg25 28_7Hz 034.mat'
    };
pivfiles = { 
    'bg25 28_7Hz 027 vortex.mat'
    'bg25 28_7Hz 028 vortex.mat'
    'bg25 28_7Hz 029 vortex.mat'
    'bg25 28_7Hz 030 vortex.mat'
    'bg25 28_7Hz 031 vortex.mat'
    'bg25 28_7Hz 032 vortex.mat'
    'bg25 28_7Hz 033 vortex.mat'
    'bg25 28_7Hz 034 vortex.mat'
    };

outfile = 'F:\Data\bg36 data/bg25 2_5BLs.csv'; % or 'F:\Data\Bg38 data/bg38 2BLs test.csv';
outmatfile = 'F:\Data\bg36 data/bg25 2_5BLs.mat'; % or 'F:\Data\Bg38 data/bg38 2BLs test.mat';

%get the noise and bias characteristics for the gyro
noise      = load_imu(noisefile);

good = noise.t >= (noise.t(end)-noise.t(1))/2;      % last half
constBiasGyro   = mean(noise.gyro(good,:),1);
beta = sqrt(3/4) * 2*max(rms(noise.gyro(good,:)));

if (length(accfiles) ~= length(kinfiles)) || (length(accfiles) ~= length(pivfiles))
    error('Different numbers of accelerometer, kinematics, or PIV files!');
end

for f = 1:length(accfiles)
    accfiles{f} = fullfile(baseaccdir,accfiles{f});
    kinfiles{f} = fullfile(basekindir,kinfiles{f});
    pivfiles{f} = fullfile(basepivdir,pivfiles{f});
end

for f = 1:length(accfiles)
    if ~exist(accfiles{f},'file')
        warning('File %s not found.  Trial will be skipped.\n', accfiles{f});
    end
    if ~exist(kinfiles{f},'file')
        warning('File %s not found.  Trial will be skipped.\n', kinfiles{f});
    end
    if ~exist(pivfiles{f},'file')
        warning('File %s not found.  Trial will be skipped.\n', pivfiles{f});
    end
    
    [~,accfile1,~] = fileparts(accfiles{f});
    [~,kinfile1,~] = fileparts(kinfiles{f});
    [~,pivfile1,~] = fileparts(pivfiles{f});
    acctok = regexp(accfile1, '[Bb]g(\d+)[ _](\d+)([a-z]?)$', 'tokens','once');
    kintok = regexp(kinfile1, '[Bb]g(\d+) [\d_]+Hz (\d+)([a-z]?)( Part \d+)?$', 'tokens','once');
    pivtok = regexp(pivfile1, '[Bb]g(\d+) [\d_]+Hz (\d+)([a-z]?) vortex$', 'tokens','once');
    
    if isempty(acctok)
        error('Cannot parse file name %s\n', accfile1);
    end
    if isempty(kintok)
        error('Cannot parse file name %s\n', kinfile1);
    end
    
    if ((str2double(acctok{1}) ~= str2double(kintok{1})) || ...
        (str2double(acctok{1}) ~= str2double(pivtok{1})))
        error('Fish numbers do not match: %s %s %s\n', accfile1, kinfile1, pivfile1);
    end
    if (str2double(acctok{2}) ~= str2double(kintok{2})) || ...
        (str2double(acctok{2}) ~= str2double(pivtok{2}))
        error('Trial numbers do not match: %s %s %s\n', accfile1, kinfile1, pivfile1);
    end
end    

nfiles = length(accfiles);

[~,acccalib] = load_imu(acccalibfile,[],'calib','axisorder',{'Y','-X','Z'},'clickzero');

load(massdistfile,'massperlen');

Kinematics = struct([]);
IMU = struct([]);
PIV = struct([]);

inputyn('','clearsaved');

if exist(outmatfile,'file') && inputyn('Continue from existing file?')
    load(outmatfile);
end
nprev = length(Kinematics);

for f = 1:nfiles
    if ~exist(kinfiles{f},'file')
        fprintf('%s not found.  Skipping.\n', kinfiles{f});
        continue;
    end
    if ~exist(accfiles{f},'file')
        fprintf('%s not found.  Skipping.\n', accfiles{f});
        continue;
    end
        
    fprintf('*** %s\n', kinfiles{f});
    if ((f > nprev) || isempty(Kinematics(f).t) || ...
            inputyn('Reload kinematics data?','default',false))
        kin1 = process_accel_kinematics(kinfiles{f},'massperlen',massperlen, ...
            'smoothdur',smoothdur, 'fishlen',fishlen);
    else
        kin1 = Kinematics(f);
    end
    if (length(kin1.t) < 20)
        warning('Very few frames (%d) digitized in file %s.  Skipping...',length(kin1.t),kinfiles{f});
        continue;
    end
    
    fprintf('*** %s\n', accfiles{f});
    if (f > nprev) || isempty(IMU(f).t) || ...
            inputyn('Reload IMU data?','default',false)
        timerange = [min(kin1.t) max(kin1.t)];
        imu1 = load_imu(accfiles{f},acccalib,'resamplerate',200, ...
            'constbiasgyro',constBiasGyro, 'timerange',timerange, ...
            'resamplerate',sampfreq);
        switch accmethod
            case 'madgwick'
                imu1 = get_orient_imu(imu1,'method','madgwick','imuposition',imuposition, ...
                    'beta',beta, 'gyroband',[0.5 10]);
        end
    else
        imu1 = IMU(f);
    end
    
    fprintf('*** %s\n', pivfiles{f});
    if ((f > nprev) || ~isfield(PIV(f),'t') || isempty(PIV(f).t) || ...
            inputyn('Reload PIV data?','default',false))
        piv1 = process_accel_piv_data(pivfiles{f}, kin1, ...
            'leftsidevortexsign',leftsidevortexsign);
    end
    
    if (f > nprev)
        fn = fieldnames(kin1);
        for i = 1:length(fn)
            Kinematics(f).(fn{i}) = kin1.(fn{i});
        end

        fn = fieldnames(imu1);
        for i = 1:length(fn)
            IMU(f).(fn{i}) = imu1.(fn{i});
        end
        
        fn = fieldnames(piv1);
        for i = 1:length(fn)
            PIV(f).(fn{i}) = piv1.(fn{i});
        end
    end
    save(outmatfile, 'Kinematics','IMU','PIV','f','-v7.3');
    
    issteady = abs(kin1.headdispfwdmn) < steadythresh;
    accstart = false(size(issteady));
    accstart(3:end) = issteady(1:end-2) & issteady(2:end-1) & ~issteady(3:end);
    accstartind = find(accstart);
    
    tailbeatnum = zeros(size(accstart));
    if ~isempty(accstartind)
        accstartind(end+1) = length(issteady)+1;
        for i = 1:length(accstartind)-1
            j = accstartind(i);
            k = 0;
            while j+k < accstartind(i+1)
                tailbeatnum(j+k) = k+1;
                k = k+1;
            end
        end
    end
    tailbeatnum(issteady) = 0;
        
    nwave = size(kin1.indpeak,2);
    nchan = 1;
    
    acc = process_accel_accel(imu1,kin1);

    [~,fn] = fileparts(kinfiles{f});
    
    out(f).filename = repmat({fn},[nchan nwave]);
    
    out(f).t = repmat(kin1.tpeak(end,:),[nchan 1]);
    out(f).tailbeatnum = repmat(tailbeatnum,[nchan 1]);
    out(f).speed = repmat(kin1.comspeedfwdmn,[nchan 1]);
    out(f).speedrms = repmat(kin1.comspeedfwdrms,[nchan 1]);
    out(f).headdisp = repmat(kin1.headdispfwdmn,[nchan 1]);
    out(f).tailamp = repmat(kin1.amp(end,:),[nchan 1]);
    out(f).freq = repmat(1./kin1.per, [nchan 1]);
    
    out(f).accfwdpk = repmat(acc.fwdpk, [nchan 1]);
    out(f).accfwdmn = repmat(acc.mean(1,:), [nchan 1]);
    out(f).accfwdiqr = repmat(acc.iqr(1,:), [nchan 1]);
    out(f).accsidemn = repmat(acc.mean(2,:), [nchan 1]);
    out(f).accsideiqr = repmat(acc.iqr(2,:), [nchan 1]);
    out(f).accupmn = repmat(acc.mean(3,:), [nchan 1]);
    out(f).accupiqr = repmat(acc.iqr(3,:), [nchan 1]);
    
    out(f).rollang = repmat(acc.orient(1,:), [nchan 1]);
    out(f).rollstd = repmat(acc.orientstd(1,:), [nchan 1]);
    out(f).pitchang = repmat(acc.orient(2,:), [nchan 1]);
    out(f).pitchstd = repmat(acc.orientstd(2,:), [nchan 1]);
    out(f).yawang = repmat(acc.orient(3,:), [nchan 1]);
    out(f).yawstd = repmat(acc.orientstd(3,:), [nchan 1]);
    
    out(f).vxcirc = repmat(piv1.vxcirc, [nchan 1]);
    out(f).vxcircstd = repmat(piv1.vxcircstd, [nchan 1]);
    out(f).vxprevcirc = repmat(piv1.vxprevcirc, [nchan 1]);
    out(f).vxprevcircstd = repmat(piv1.vxprevcircstd, [nchan 1]);
    out(f).vxdist = repmat(piv1.vxdist, [nchan 1]);
    out(f).vxdiststd = repmat(piv1.vxdiststd, [nchan 1]);
end

putvar out;
save_struct_as_table(outfile, out);
