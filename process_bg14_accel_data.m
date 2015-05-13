function process_bg14_accel_data

accmethod = 'old';
smoothdur = 0.5;
doplot = true;
dodiagnostic = true;

imuposition = [6.6 11.4 -7];
emgposition.L1 = [1 0.42];
emgposition.L2 = [1 0.55];
emgposition.L3 = [1 0.68];
emgposition.L4 = [1 0.77];
emgposition.R1 = [2 0.42];
emgposition.R2 = [2 0.55];
emgposition.R3 = [2 0.68];
emgposition.R4 = [2 0.77];

threshold = [...
   -0.0693   -0.1221   -0.2544   -0.1476   -0.2329   -0.2227   -0.0895   -0.1222; ...
    0.0652    0.1203    0.2070    0.1311    0.2462    0.2650    0.1102    0.1441];
interburstdur = [0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];
minspikes = [1 1 3 3 1 3 1 2];

steadythresh = 0.03;

acccalibfile = 'rawdata/bg14/Accelerometer/calib001.h5';
massdistfile = 'rawdata/fishmass.mat';

emgfiles = {'rawdata/bg14/EMG/bg14_005.mat'};
accfiles = {'rawdata/bg14/Accelerometer/bg14_005.h5'};
kinfiles = {'rawdata/bg14/digitizeFish/Bg14_005.mat'};
outfile = 'rawdata/bg14/bg14data.csv';

nfiles = length(emgfiles);

[~,acccalib] = load_imu(acccalibfile,[],'calib','axisorder',{'Y','-Z','-X'});

load(massdistfile,'massperlen');

for f = 1:nfiles
    imu = load_imu(accfiles{f},acccalib,'resamplerate',200);
    switch accmethod
        case 'old'
            imu = get_orient_imu_updated(imu,'imuposition',imuposition, ...
                'getoffset',false,'gyrooffset',[-30 -29]);
    end
    
    emg = importLabChart(emgfiles{f},[],'outformat','new');
    
    kin = process_accel_kinematics(kinfiles{f},'massperlen',massperlen, ...
        'smoothdur',smoothdur);
    tend = kin.t(end);
    kin.t = kin.t - tend;
    kin.tpeak = kin.tpeak - tend;
    kin.tcurvepeak = kin.tcurvepeak - tend;

    issteady = abs(kin.headdispfwdmn) < steadythresh;
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
    
    emg = process_accel_emg(emg,kin, 'spikethreshold',threshold, ...
        'interburstdur',interburstdur, 'minspikes',minspikes, 'goodchan',true(size(minspikes)), ...
        'emgposition',emgposition);
    
    emgpos1 = struct2cell(emgposition);
    emgpos1 = cat(1,emgpos1{:});
    nwave = size(kin.tpeak,2);
    nchan = size(emgpos1,1);
    
    emgpos = repmat(emgpos1(:,2),[1 nwave]);
    emgside = repmat(emgpos1(:,1),[1 nwave]);

    if dodiagnostic
        figureseries('Burst data');
        clf;
        
        subplot(2,1,1);
        h = plotgroups(emg.burstctr, emg.burstdur, {emgside,emgpos}, {'f','cm'});
        axis tight;
        ylabel('Burst duration (sec)');
        legend(h,fieldnames(emgposition),'location','best');
        
        bd = emg.burstdur;
        medburstdur = zeros(1,nwave);
        k = 3:nwave-2;
        medburstdur(3:end-2) = nanmedian(cat(1,bd(:,k-2),bd(:,k-1),bd(:,k),bd(:,k+1),bd(:,k+2)));
        medburstdur(1) = nanmedian(cat(1,bd(:,1),bd(:,2),bd(:,3)));
        medburstdur(2) = nanmedian(cat(1,bd(:,1),bd(:,2),bd(:,3),bd(:,4)));
        medburstdur(end-1) = nanmedian(cat(1,bd(:,end-3),bd(:,end-2),bd(:,end-1),bd(:,end)));
        medburstdur(end) = nanmedian(cat(1,bd(:,end-2),bd(:,end-1),bd(:,end)));
        
        medburstdur = repmat(medburstdur,[nchan 1]);
        %CONTINUE here
        
        isoutlier = abs(bd - medburstdur)./medburstdur > 0.5;
        addplot(emg.burstctr(isoutlier), emg.burstdur(isoutlier), 'ro', ...
            'MarkerSize',12);
        
        subplot(2,1,2);
        plotgroups(emg.burstctr, emg.burstfreq, {emgside,emgpos}, {'f','cm'});
        axis tight;
        ylabel('Burst frequency (Hz)');
        xlabel('Time (sec)');
    end
    
    acc = process_accel_accel(imu,kin);


    if doplot
        accfwdpk = acc.fwdpk;
        accfwdpk = repmat(accfwdpk,[size(emgpos,1) 1]);
        
        [~,ord] = sort(acc.fwdpk);
        jit1 = zeros(size(acc.fwdpk));
        %this gives us the rank of each acceleration
        jit1(ord) = 1:length(acc.fwdpk);
        jit1 = (jit1-1)/range(jit1) * 0.02;
        jit1 = repmat(jit1,[size(emgpos,1) 1]);
        
        ptsz = emg.burstamp;
        ptsz = (ptsz - min(ptsz(:))) / range(ptsz(:)) * 12 + 4;
        ptsz = pi*ptsz.^2;
        
        isleft = emgside == 1;
        X = emgpos(isleft)+jit1(isleft)-0.011;
        
        clf;
        hold on;
        plot([X(:) X(:)]', [emg.burstonphase(isleft) emg.burstoffphase(isleft)]','k-');        
        scatter(X(:), emg.burstonphase(isleft), ptsz(isleft), accfwdpk(isleft), 'filled','o');
        scatter(X(:), emg.burstoffphase(isleft), ptsz(isleft), accfwdpk(isleft), 'filled','o');
        
        X = emgpos(~isleft)+jit1(~isleft)+0.011;
        plot([X(:) X(:)]', [emg.burstonphase(~isleft) emg.burstoffphase(~isleft)]', 'r-');
        scatter(X(:), emg.burstonphase(~isleft), ptsz(~isleft), accfwdpk(~isleft), 'filled','s');
        scatter(X(:), emg.burstoffphase(~isleft), ptsz(~isleft), accfwdpk(~isleft), 'filled','s');
        hold off;
    end
    
    emgsidec = repmat('L',size(emgside));
    emgsidec(emgside == 2) = 'R';
    
    [~,fn] = fileparts(accfiles{f});
    
    out.filename = repmat({fn},[nchan nwave]);
    
    out.emgpos = emgpos;
    out.emgside = emgsidec;
    out.t = repmat(kin.tpeak(end,:),[nchan 1]);
    out.tailbeatnum = repmat(tailbeatnum,[nchan 1]);
    out.speed = repmat(kin.comspeedfwdmn,[nchan 1]);
    out.speedrms = repmat(kin.comspeedfwdrms,[nchan 1]);
    out.headdisp = repmat(kin.headdispfwdmn,[nchan 1]);
    out.headamp = repmat(kin.amp(1,:),[nchan 1]);
    out.tailamp = repmat(kin.amp(1,:),[nchan 1]);
    out.freq = repmat(1./kin.per, [nchan 1]);
    out.wavespeed = repmat(kin.wavespeed, [nchan 1]);
    out.wavelen = repmat(nanmean(kin.wavelen), [nchan 1]);
    
    out.accfwdpk = repmat(acc.fwdpk, [nchan 1]);
    out.accfwdmn = repmat(acc.mean(1,:), [nchan 1]);
    out.accfwdiqr = repmat(acc.iqr(1,:), [nchan 1]);
    out.accsidemn = repmat(acc.mean(2,:), [nchan 1]);
    out.accsideiqr = repmat(acc.iqr(2,:), [nchan 1]);
    out.accupmn = repmat(acc.mean(3,:), [nchan 1]);
    out.accupiqr = repmat(acc.iqr(3,:), [nchan 1]);
    
    out.rollang = repmat(acc.orient(1,:), [nchan 1]);
    out.rollstd = repmat(acc.orientstd(1,:), [nchan 1]);
    out.pitchang = repmat(acc.orient(2,:), [nchan 1]);
    out.pitchstd = repmat(acc.orientstd(2,:), [nchan 1]);
    out.yawang = repmat(acc.orient(3,:), [nchan 1]);
    out.yawstd = repmat(acc.orientstd(3,:), [nchan 1]);
    
    out.burston = emg.burston;
    out.burstoff = emg.burstoff;
    out.burstonphase = emg.burstonphase;
    out.burstoffphase = emg.burstoffphase;
    out.burstduty = emg.burstduty;
    out.burstamp = emg.burstamp;
    out.burstint = emg.burstint;
    
    save_struct_as_table(outfile, out);
end