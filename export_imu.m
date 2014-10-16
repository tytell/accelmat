function export_imu(filename,calibfile, outfile, varargin)

opt.resamplerate = [];
opt.kinematicsfile = '';
opt.emgfile = '';
opt.outputformat = '';
opt = parsevarargin(opt,varargin, 4);

if nargin == 0
    [fn,pn] = uigetfile('*.h5','Choose IMU data file');
    filename = fullfile(pn,fn);
    
    [fn,pn] = uigetfile('*.mat','Choose calibration file',pn);
    calibfile = fullfile(pn,fn);
    
    [fn,pn1] = uigetfile('*.mat','Choose kinematics file (or Cancel for none)',pn);
    if ischar(fn) && ~isempty(fn)
        kinfile = fullfile(pn1,fn);
    else
        kinfile = '';
    end
    
    [fn,pn1] = uigetfile('*.mat','Choose EMG file (or Cancel for none)',pn);
    if ischar(fn) && ~isempty(fn)
        emgfile = fullfile(pn1,fn);
    else
        emgfile = '';
    end
    
    [fn,pn] = uiputfile('*.txt','Choose output file',pn);
    outfile = fullfile(pn,fn);
else
    kinfile = opt.kinematicsfile;
    emgfile = opt.emgfile;
end;

if (ischar(calibfile) && exist(calibfile,'file'))
    load(calibfile,'calib');
elseif isstruct(calibfile)
    calib = calibfile;
end
imu = load_imu(filename,calib,'resample',false);

imurate = 1./(imu.t(2) - imu.t(1));
resampletxt = sprintf('IMU=%.2fHz', imurate);
defaultrate = imurate;

if ~isempty(kinfile)
    K = load(kinfile,'haxmmss','haymmss','humms','hvmms','hxmm',...
        'hymm','txmm','tymm','t');
    kinrate = 1/(K.t(2) - K.t(1));
    resampletxt = [resampletxt sprintf(', Kinematics=%dHz', kinrate)];
end

if ~isempty(emgfile)
    E = importLabChart(emgfile,'');
    if ~isfield(E,'t')
        error('Cannot handle multirate LabChart files yet...');
    end
    emgrate = E.samplerate(1);
    resampletxt = [resampletxt sprintf(', EMG=%dHz', emgrate)];
    defaultrate = emgrate;
end

if isempty(opt.resamplerate)
    fprintf('Resample rate?\n  (%s) default=%d', resampletxt,defaultrate);
    r = input(': ');
    if ~isempty(r) && isnumeric(r)
        resamplerate = r;
    else
        resamplerate = defaultrate;
        fprintf('Resampling at %.2fHz.\n', defaultrate);
    end
else
    resamplerate = opt.resamplerate;
end

imu = get_orient_imu(imu, 'getoffset');

t0 = imu.t;
t1 = (t0(1):1/resamplerate:t0(end))';
X0 = [imu.orient imu.acchi];
X0(:,1:3) = X0(:,1:3) * 180/pi;     %convert to degrees
X0(:,4:6) = X0(:,4:6) * 9.800;        %convert to m/s^2
X1 = zeros(length(t1),6);
for i = 1:size(X0,2)
    X1(:,i) = interp1(t0,X0(:,i), t1);
end
X1 = [t1 X1];

chantitle = {'Pitch','Roll','Yaw','UpAcc','FwdAcc','SideAcc'};
tplt = '%.5f\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.5f';
if ~isempty(kinfile)
    X2 = [K.hxmm' K.hymm' K.txmm' K.tymm' K.humms' K.hvmms' K.haxmmss' K.haymmss'];
    K.t = K.t - K.t(end);
    Xkin = zeros(length(t1),8);
    for i = 1:size(X2,2)
        Xkin(:,i) = interp1(K.t,X2(:,i), t1);
    end
    Xkin = Xkin / 1000;       % convert to m
    
    X1 = [X1 Xkin];
    
    chantitle = [chantitle, {'FwdPos','SidePos','TailFwd','TailSide','FwdVel','SideVel',...
        'FwdAccVid','SideAccVid'}];
    tplt = [tplt '\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f'];
end

annote.event = {};
annote.starttime = [];
annote.duration = [];

if ~isempty(emgfile)
    emgchantitle = fieldnames(E);
    trigchanind = find(ismember(lower(emgchantitle),{'trigger','trig'}));
    if ~isempty(trigchanind)
        trigchan = emgchantitle{trigchanind};
        t0 = first(E.t, E.(trigchan) < 3);
    else
        t0 = E.t(end);
    end
    
    emgt = E.t - t0;
    
    if isfield(E,'commenttxt')
        annote.event = E.commenttxt;
        annote.starttime = E.commentt;
        annote.duration = zeros(size(E.commentt));
    end

    len0 = length(E.t);
    E = struct2cell(E);
    len = cellfun(@length,E);
    
    ischan = (len == len0) & ~ismember(emgchantitle,{'t'});
    if ~isempty(trigchanind)
        ischan(trigchanind) = false;
    end
    X2 = cat(2,E{ischan});
    
    Xemg = zeros(length(t1),sum(ischan));
    for i = 1:size(X2,2)
        Xemg(:,i) = interp1(emgt,X2(:,i), t1);
    end
    
    X1 = [X1 Xemg];
    
    emgchantitle = emgchantitle(ischan);
    chantitle = [chantitle emgchantitle(:)'];
    tplt = [tplt repmat('\t%.5f',[1 sum(ischan)])];
    
    
end

if isempty(opt.outputformat)
    fprintf('Output file format:\n');
    fprintf(' T. LabChart Text\n M. Matlab\n A. All (default)\n');
    fmt1 = input('Choose format: ','s');
    
    switch lower(fmt1)
        case {'t','text'}
            fmt = {'text'};
        case {'m','mat','matlab'}
            fmt = {'mat'};
        case {'a','all'}
            fmt = {'text','mat'};
        otherwise
            fmt = {'text','mat'};
    end
else
    fmt = lower(opt.outputformat);
    if ~iscell(fmt)
        fmt = {fmt};
    end
    if any(~ismember(fmt,{'text','mat','edf'}))
        error('Unrecognized output format');
    end
end

if ismember('text',fmt)
    fid = fopen(outfile,'w');
    titletplt = ['ChannelTitle=\t' repmat('%s\t',[1 length(chantitle)-1]) '%s\n'];
    fprintf(fid, titletplt, chantitle{:});
    fprintf(fid, [tplt '\n'], X1');
    fclose(fid);
end

[pn,fn,ext] = fileparts(outfile);

if ismember('mat',fmt)
    outmatfile = fullfile(pn,[fn '.mat']);
    
    C = mat2cell(X1, length(t1), ones(1,size(X1,2)));
    S = cell2struct(C(:),['t' chantitle],1);
    save(outmatfile,'-struct','S');
end

if ismember('edf',fmt)
    outedffile = fullfile(pn,[fn '.edf']);
    
    header.samplerate = resamplerate;
    header.duration = length(t1)/resamplerate;
    header.labels = chantitle;
    header.annotation = annote;
    SaveEDF(outedffile,X1(:,2:end),header);
end

    
            
        




    