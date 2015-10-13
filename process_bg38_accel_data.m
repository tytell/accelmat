function process_bg38_accel_data

accmethod = 'madgwick';
smoothdur = 0.5;
leftsidevortexsign = -1;

doplot = true;
dodiagnostic = true;

fishlen = 130;       % mm

sampfreq = 100;

imuposition = [6.6 11.4 -7];

steadythresh = 0.03;

acccalibfile = 'rawdata/bg38/calib001.h5';
massdistfile = './fishmass.mat';

accfiles = {
    'rawdata/bg38/bg38_011.h5'
    'rawdata/bg38/bg38_015.h5'
    'rawdata/bg38/bg38_019.h5'
    'rawdata/bg38/bg38_023.h5'
    };
kinfiles = { 
    'rawdata/bg38/bg38 16_5Hz 011.mat'
    'rawdata/bg38/bg38 16_5Hz 015.mat'
    'rawdata/bg38/bg38 21_3Hz 019.mat'
    'rawdata/bg38/bg38 21_3Hz 023.mat'
    };
pivfiles = { 
    'rawdata/bg38/bg38 16_5Hz 011 vortex.mat'
    'rawdata/bg38/bg38 16_5Hz 015 vortex.mat'
    'rawdata/bg38/bg38 21_3Hz 019 vortex.mat'
    'rawdata/bg38/bg38 21_3Hz 023 vortex.mat'
    };

outfile = 'rawdata/bg38/bg38data.csv';
outmatfile = 'rawdata/bg38/bg38data.mat';

%get the noise and bias characteristics for the gyro
noise      = load_imu('rawdata/bg38/noisytest001.h5');

good = noise.t >= (noise.t(end)-noise.t(1))/2;      % last half
constBiasGyro   = mean(noise.gyro(good,:),1);
beta = sqrt(3/4) * 2*max(rms(noise.gyro(good,:)));

for f = 1:min([length(accfiles) length(kinfiles)])
    if ~exist(accfiles{f},'file')
        warning('File %s not found.  Trial will be skipped.\n', accfiles{f});
    end
    if ~exist(kinfiles{f},'file')
        warning('File %s not found.  Trial will be skipped.\n', kinfiles{f});
    end
    
    [~,accfile1,~] = fileparts(accfiles{f});
    [~,kinfile1,~] = fileparts(kinfiles{f});
    acctok = regexp(accfile1, '[Bb]g(\d+)[ _](\d+)([a-z]?)$', 'tokens','once');
    kintok = regexp(kinfile1, '[Bb]g(\d+)[ _].*[ _\-](\d+)([a-z]?)$', 'tokens','once');
    
    if isempty(acctok)
        error('Cannot parse file name %s\n', accfile1);
    end
    if isempty(kintok)
        error('Cannot parse file name %s\n', kinfile1);
    end
    
    if (str2double(acctok{1}) ~= str2double(kintok{1}))
        error('Fish numbers do not match: %s %s\n', accfile1, kinfile1);
    end
    if (str2double(acctok{2}) ~= str2double(kintok{2}))
        error('Trial numbers do not match: %s %s\n', accfile1, kinfile1);
    end
end    

nfiles = length(accfiles);

[~,acccalib] = load_imu(acccalibfile,[],'calib','axisorder',{'Y','Z','-X'});

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
    
    if ((f > nprev) || isempty(PIV(f).t) || ...
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
    out(f).vxprevcirc = repmat(piv1.vxcirc, [nchan 1]);
    out(f).vxprevcircstd = repmat(piv1.vxcircstd, [nchan 1]);
    out(f).vxdist = repmat(piv1.vxcirc, [nchan 1]);
    out(f).vxdiststd = repmat(piv1.vxcircstd, [nchan 1]);
end

putvar out;
save_struct_as_table(outfile, out);
