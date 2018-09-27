function process_bg14_accel_data

accmethod = 'madgwick';   % shouldn't need changing
smoothdur = 0.5;% shouldn't need changing
doplot = true;          % shouldn't need changing
dodiagnostic = true;    % shouldn't need changing

% shouldn't need changing
nrunmedian = 7;
duroutlierfrac = 0.8;
freqoutlierfrac = 0.2;

% update per fish
imuposition = [11.4 -5.0 13.9];        % relative to COM, x, y, z coord of IMU in mm
emgposition.L1 = [1 0.349];      %white muscle     % 1 = L, 2 = R, and then position in fraction of L
emgposition.L2 = [1 0.566];
emgposition.L3 = [1 0.612];    %pulled out after 007
emgposition.L4 = [1 0.717];    %pulled out after 007            
emgposition.R1 = [2 0.434]; 
emgposition.R2 = [2 0.546];
emgposition.R3 = [2 0.658];
emgposition.R4 = [2 0.763];        


% defaults for burst detection - shouldn't need to change
threshold = [...
   -0.0693   -0.1221   -0.2544   -0.1476   -0.2329   -0.2227   -0.0895   -0.1222; ...
    0.0652    0.1203    0.2070    0.1311    0.2462    0.2650    0.1102    0.1441];
interburstdur = [0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];
minspikes = [1 1 3 3 1 3 1 2];
goodchan = true(size(minspikes));

% fraction of body length movement for beginning of acceleration
steadythresh = 0.03;

acccalibfile = 'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/calib001.h5';
noisefile = 'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/noisytest001.h5';

massdistfile = 'fishmass.mat';

emgfiles = {'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_001.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_002.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_003.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_003.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_004.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_005.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_006.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_007.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_008.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_009.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_010.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_011.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_011.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_011.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_012.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_013.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_014.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_014.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_014.mat'...     % 19
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_014.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_015.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_015.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_016.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_016.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_016.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_016.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_016.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_017.mat'...     % 28
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_017.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_018.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_018.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_019.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_020.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_021.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_022.mat'...     % 35
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_023.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_024.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_025.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_026.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_027.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_027.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_028.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_028.mat'...     % 43
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_029.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_031.mat'...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/EMGs/bg14_032.mat'};
accfiles = {'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_001.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_002.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_003.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_003.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_004.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_005.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_006.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_007.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_008.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_009.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_010.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_011.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_011.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_011.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_012.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_013.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_014.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_014.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_014.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_014.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_015.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_015.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_016.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_016.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_016.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_016.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_016.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_017.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_017.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_018.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_018.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_019.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_020.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_021.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_022.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_023.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_024.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_025.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_026.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_027.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_027.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_028.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_028.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_029.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_031.h5' ...
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Accelerometer/bg14_032.h5'};
kinfiles = {'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_001.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_002.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_003a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_003b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_004.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_005.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_006.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_007.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_008.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_009.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_010.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_011a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_011b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_011c.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_012.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_013.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_014a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_014b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_014c.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_014d.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_015a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_015b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_016a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_016b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_016c.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_016d.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_016e.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_017a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_017b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_018a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_018b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_019.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_020.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_021.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_022.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_023.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_024.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_025.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_026.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_027a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_027b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_028a.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_028b.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_029.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_031.mat'
    'T:/Bluegill Acceleration/2015-03-26 Bg14/Bg14 digitizeFish files/Bg14_032.mat'};

outfile = 'T:/Bluegill Acceleration/2015-03-26 Bg14/bg14data.csv';
outmatfile = 'T:/Bluegill Acceleration/2015-03-26 Bg14/bg14data.mat';

%don't change after here
nfiles = length(emgfiles);

[~,acccalib] = load_imu(acccalibfile,[],'calib','clickzero'); %'axisorder',{'-Y','-Z','X'});

%get the noise and bias characteristics for the gyro
noise      = load_imu(noisefile);

good = noise.t >= (noise.t(end)-noise.t(1))/2;      % last half
constBiasGyro   = mean(noise.gyro(good,:),1);
beta = sqrt(3/4) * 2*max(rms(noise.gyro(good,:)));

load(massdistfile,'massperlen');

Kinematics = struct([]);
EMG = struct([]);
IMU = struct([]);

inputyn('','clearsaved');

if exist(outmatfile,'file') && inputyn('Continue from existing file?')
    load(outmatfile);
end
nprev = length(Kinematics);

for f = 1:nfiles
    fprintf('*** %s\n', kinfiles{f});
    if ((f > nprev) || isempty(Kinematics(f).t) || ...
            inputyn('Reload kinematics data?','default',false))
        kin1 = process_accel_kinematics(kinfiles{f},'massperlen',massperlen, ...
            'smoothdur',smoothdur);
    else
        kin1 = Kinematics(f);
    end
    
    if (f > nprev) || isempty(IMU(f).t) || ...
            inputyn('Reload IMU data?','default',false)
        timerange = [min(kin1.t) max(kin1.t)];
        imu1 = load_imu(accfiles{f},acccalib,'resamplerate',200, ...
            'constbiasgyro',constBiasGyro, 'timerange',timerange);
        switch accmethod
            case 'madgwick'
                imu1 = get_orient_imu(imu1,'method','madgwick','imuposition',imuposition, ...
                    'beta',beta, 'gyroband',[0.5 10]);
        end
    else
        imu1 = IMU(f);
    end
    
    if (f > nprev) || isempty(EMG(f).t) || ...
            inputyn('Process EMG data again?','default',false)

        emg1 = importLabChart(emgfiles{f},[],'outformat','new');
        
        if (f <= nprev) && ~isempty(EMG(f).t)
            override = EMG(f).burstoverride;
            threshold = EMG(f).spikethreshold;
            interburstdur = EMG(f).interburstdur;
            minspikes = EMG(f).minspikes;
            goodchan = EMG(f).goodchan;            
        else
            override = struct([]);
        end
        emg1 = process_accel_emg(emg1,kin1, 'spikethreshold',threshold, ...
            'interburstdur',interburstdur, 'minspikes',minspikes, 'goodchan',goodchan, ...
            'emgposition',emgposition,'override',override);
        
        %update for the next run
        threshold = emg1.spikethreshold;
        interburstdur = emg1.interburstdur;
        minspikes = emg1.minspikes;
        goodchan = emg1.goodchan;
    else
        emg1 = EMG(f);
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
        
        fn = fieldnames(emg1);
        for i = 1:length(fn)
            EMG(f).(fn{i}) = emg1.(fn{i});
        end
    end
    save(outmatfile, 'Kinematics','EMG','IMU','f','-v7.3');
    
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
    nchan = size(emg1.pos,1);
    
    emgpos = repmat(emg1.pos,[1 nwave]);
    emgside = repmat(emg1.side,[1 nwave]);

    if dodiagnostic
        figureseries('Burst data');
        clf;
        
        hax = zeros(2,1);
        hax(1) = subplot(2,1,1);
        h = plotgroups(emg1.burstctr, emg1.burstdur, {emgside,emgpos}, {'f','cm'});
        axis tight;
        ylabel('Burst duration (sec)');
        legend(h,fieldnames(emgposition),'location','best');
        
        nmed = floor(nrunmedian/2);
        
        medburstdur = NaN(2*nmed+1,nchan,nwave);
        k = 1:nwave;
        off = -nmed:nmed;
        for i = 1:length(off)
            koff = k+off(i);
            good = (koff >= 1) & (koff <= nwave);
            medburstdur(i,:,k(good)) = emg1.burstdur(:,koff(good));
        end
        medburstdur = flatten(medburstdur,1:2);
        medburstdur = nanmedian(medburstdur);
        medburstdur = repmat(medburstdur,[nchan 1]);
        
        isoutlier = abs(emg1.burstdur - medburstdur)./medburstdur > duroutlierfrac;
        addplot(emg1.burstctr(isoutlier), emg1.burstdur(isoutlier), 'ro', ...
            'MarkerSize',12);
        
        hax(2) = subplot(2,1,2);
        plotgroups(emg1.burstctr, emg1.burstfreq, {emgside,emgpos}, {'f','cm'});
        axis tight;
        ylabel('Burst frequency (Hz)');
        xlabel('Time (sec)');
        
        medburstfreq = NaN(2*nmed+1,nchan,nwave);
        k = 1:nwave;
        off = -nmed:nmed;
        for i = 1:length(off)
            koff = k+off(i);
            good = (koff >= 1) & (koff <= nwave);
            medburstfreq(i,:,k(good)) = emg1.burstfreq(:,koff(good));
        end
        medburstfreq = flatten(medburstfreq,1:2);
        medburstfreq = nanmedian(medburstfreq);
        medburstfreq = repmat(medburstfreq,[nchan 1]);
        
        isoutlier = abs(emg1.burstfreq - medburstfreq)./medburstfreq > freqoutlierfrac;
        addplot(emg1.burstctr(isoutlier), emg1.burstfreq(isoutlier), 'ro', ...
            'MarkerSize',12);
        
        linkaxes(hax,'x');
    end
    
    acc = process_accel_accel(imu1,kin1);

    if doplot
        jit1 = linspace(0,1, size(emg1.burstonphase,2)) * 0.02;
        jit1 = repmat(jit1,[nchan 1]);
        
        isleft = emgside == 1;
        X = emgpos(isleft)+jit1(isleft)-0.011;
        addplot([X(:) X(:)]', [emg1.burstonphase(isleft) emg1.burstoffphase(isleft)]', 'bo-');
        
        X = emgpos(~isleft)+jit1(~isleft)+0.011;
        addplot([X(:) X(:)]', [emg1.burstonphase(~isleft) emg1.burstoffphase(~isleft)]', 'gs-');
    
    end
    
    emgsidec = repmat('L',size(emgside));
    emgsidec(emgside == 2) = 'R';
    
    [~,fn] = fileparts(kinfiles{f});
    
    out(f).filename = repmat({fn},[nchan nwave]);
    
    out(f).emgpos = emgpos;
    out(f).emgside = emgsidec;
    out(f).t = repmat(kin1.tpeak(end,:),[nchan 1]);
    out(f).tailbeatnum = repmat(tailbeatnum,[nchan 1]);
    out(f).speed = repmat(kin1.comspeedfwdmn,[nchan 1]);
    out(f).speedrms = repmat(kin1.comspeedfwdrms,[nchan 1]);
    out(f).headdisp = repmat(kin1.headdispfwdmn,[nchan 1]);
    out(f).headamp = repmat(kin1.amp(1,:),[nchan 1]);
    out(f).tailamp = repmat(kin1.amp(end,:),[nchan 1]);
    out(f).freq = repmat(1./kin1.per, [nchan 1]);
    out(f).wavespeed = repmat(kin1.wavespeed, [nchan 1]);
    out(f).wavelen = repmat(nanmean(kin1.wavelen), [nchan 1]);
    
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
    
    out(f).burston = emg1.burston;
    out(f).burstoff = emg1.burstoff;
    out(f).burstonphase = emg1.burstonphase;
    out(f).burstoffphase = emg1.burstoffphase;
    out(f).burstduty = emg1.burstduty;
    out(f).burstamp = emg1.burstamp;
    out(f).burstint = emg1.burstint;
end

save_struct_as_table(outfile, out);
