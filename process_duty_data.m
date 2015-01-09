function process_duty_data

naccelbeats = 4;

if (~getvar('duty','accpk','tailamp','per','cycletype','speed',...
    'dutyerr','accpkerr','tailamperr','pererr') || ...
    inputyn('Redo calculations?','default',false))

    files = {    'rawdata/Bg9 Export_imu/Bg9_003.mat',...
        'rawdata/Bg9 Export_imu/Bg9_004.mat',...
        'rawdata/Bg9 Export_imu/Bg9_005.mat',...
        'rawdata/Bg9 Export_imu/Bg9_006.mat',...
        'rawdata/Bg9 Export_imu/Bg9_007.mat',...
        'rawdata/Bg9 Export_imu/Bg9_008.mat',...
        'rawdata/Bg9 Export_imu/Bg9_009.mat',...
        'rawdata/Bg9 Export_imu/Bg9_010.mat',...
        'rawdata/Bg9 Export_imu/Bg9_011.mat',...
        'rawdata/Bg9 Export_imu/Bg9_012.mat',...
        'rawdata/Bg9 Export_imu/Bg9_013.mat',...
        'rawdata/Bg9 Export_imu/Bg9_014.mat',...
        'rawdata/Bg9 Export_imu/Bg9_015.mat',...
        'rawdata/Bg9 Export_imu/Bg9_016.mat',...
        'rawdata/Bg9 Export_imu/Bg9_018.mat',...
        'rawdata/Bg9 Export_imu/Bg9_019.mat',...
        'rawdata/Bg9 Export_imu/Bg9_020.mat',...
        'rawdata/Bg9 Export_imu/Bg9_021.mat',...
        'rawdata/Bg9 Export_imu/Bg9_022.mat',...
        'rawdata/Bg9 Export_imu/Bg9_024.mat',...
        'rawdata/Bg9 Export_imu/Bg9_025.mat',...
        'rawdata/Bg9 Export_imu/Bg9_027.mat',...
        'rawdata/Bg9 Export_imu/Bg9_028.mat',...
        'rawdata/Bg9 Export_imu/Bg9_029.mat'};

    datafile = 'rawdata/Bg9 Export_imu/datafiles.mat';
    burstdetectfile = 'burstdetect.mat';

    load(datafile);

    isacc = regexp(content,'[Aa]ccel', 'match','once');
    isacc = cellfun(@(x) ~isempty(x), isacc);

    issteady = regexp(content,'[Ss]teady', 'match','once');
    issteady = cellfun(@(x) ~isempty(x), issteady);

    ismixed = isacc & issteady;

    isacc(ismixed) = false;
    issteady(ismixed) = false;

    duty = [];
    dutyerr = [];
    cycletype = [];
    accpk = [];
    accpkerr = [];
    tailamp = [];
    tailamperr = [];
    per = [];
    pererr = [];
    speed2 = [];
    for i = 1:length(files)
        fprintf('%s...\n', files{i});
        [pn,fn] = fileparts(files{i});
        numstr = regexp(fn,'\d{3,3}','match');

        filenum1 = str2double(numstr);

        S1 = process_accel_data(files{i}, 'burstparams',burstdetectfile);
        nchan = size(S1.burstduty,2);

        k = find(filenum == filenum1);
        for j = 1:length(k)
            t1 = tstart(k(j));
            t2 = tend(k(j));

            ind = find((S1.tbeat >= t1) & (S1.tbeat <= t2));
            if isacc(k(j))
                duty1 = NaN(naccelbeats,nchan);
                accpk1 = NaN(naccelbeats,1);
                tailamp1 = NaN(naccelbeats,1);
                per1 = NaN(naccelbeats,1);
                cycletype1 = NaN(naccelbeats,1);

                if (length(ind) > naccelbeats)
                    ind = ind(1:naccelbeats);
                end
                good = false(naccelbeats,1);
                good(1:length(ind)) = true;

                duty1(good,:) = S1.burstduty(ind,:);
                dutyerr1 = NaN(size(duty1));
                accpk1(good,:) = S1.accpk(ind,:);
                accpkerr1 = NaN(size(accpk1));
                tailamp1(good,:) = S1.tailamp(ind,:);
                tailamperr1 = NaN(size(tailamp1));
                per1(good,:) = S1.per(ind,:);
                pererr1 = NaN(size(per1));
                cycletype1(good) = 1:length(ind);
            elseif issteady(k(j))
                duty1 = nanmedian(S1.burstduty(ind,:));
                dutyerr1 = iqr(S1.burstduty(ind,:));
                accpk1 = nanmedian(S1.accpk(ind,:));
                accpkerr1 = iqr(S1.accpk(ind,:));
                tailamp1 = nanmedian(S1.tailamp(ind,:));
                tailamperr1 = iqr(S1.tailamp(ind,:));
                per1 = nanmedian(S1.per(ind,:));
                pererr1 = iqr(S1.per(ind,:));
                cycletype1 = 0;
            end 
            duty = [duty; duty1];
            dutyerr = [dutyerr; dutyerr1];
            accpk = [accpk; accpk1];
            accpkerr = [accpkerr; accpkerr1];
            tailamp = [tailamp; tailamp1];
            tailamperr = [tailamperr; tailamperr1];
            per = [per; per1];
            pererr = [pererr; pererr1];
            cycletype = [cycletype; cycletype1];
            speed2 = [speed2; speed(k(j)) * ones(size(cycletype1))];
        end
    end

    speed = speed2;
    putvar duty accpk tailamp per cycletype speed ...
        dutyerr accpkerr tailamperr pererr;
end

chanpos = [1 2 2 3 3];
chanside = 'RLRLR';

cmap1 = hsv2rgb([0 0 0; 0.5 0.9 0.7; 0.8 0.9 0.7]);
colormap(cmap1);

figureseries('Duty1');
h = plotgroups(repmat(cycletype,[1 5]),duty,...
    {repmat(chanpos,[size(cycletype,1) 1])},{'cmf'},...
    'means','traces','on','error','std','xoff',0.08, ...
    'MarkerSize',20,'errLineWidth',3);
legend(h,'anterior','posterior','peduncle');
xtick(0:4,{'steady','1','2','3','4'});
set(gca,'TickDir','out');
xlabel('Tail beat');
ylabel('Duty cycle');

figureseries('Duty2');
colormap(cmap1);
h = plotgroups(repmat(cycletype>0,[1 5]),duty,...
    {repmat(chanpos,[size(cycletype,1) 1])},{'cmf'},...
    'means','traces','on','error','std','xoff',0.08, ...
    'MarkerSize',20,'errLineWidth',3);
legend(h,'anterior','posterior','peduncle');
xtick(0:1,{'steady','accel'});
set(gca,'TickDir','out');
ylabel('Duty cycle');

figureseries('Duty3');
cmap2 = hsv2rgb([0.6 0.9 0.7; 0.3 0.9 0.7; 0 0.9 0.7]);
colormap(cmap2);
good = speed < 2;
h = plotgroups(repmat(cycletype(good)>0,[1 2]),duty(good,4:5),...
    {repmat(speed(good),[1 2])},{'cmf'},...
    'means','traces','on','error','std','xoff',0.08, ...
    'MarkerSize',20,'errLineWidth',3,'legend');
xtick(0:1,{'steady','accel'});
ylabel('Duty cycle');
set(gca,'TickDir','out');

figureseries('Accel');
plotgroups(cycletype,accpk/1000,'means','error','std',...
    'color','k','markers','s','MarkerSize',20,'errLineWidth',3);
xtick(0:4,{'steady','1','2','3','4'});
set(gca,'TickDir','out');
xlabel('Tail beat');
ylabel('Acceleration (m s-1)');
xlim([-0.5 4.5]);

anovan(duty(:),{flatten(repmat(cycletype>0,[1 5])) flatten(repmat(speed,[1 5])) ...
    flatten(repmat(chanpos,[size(cycletype,1) 1]))},'model','full',...
    'varnames',{'type','speed','pos'});

good = speed < 2;
anovan(flatten(duty(good,4:5)),{flatten(repmat(cycletype(good)>0,[1 2])) ...
    flatten(repmat(speed(good),[1 2]))},...
    'model','full',...
    'varnames',{'type','speed'});

