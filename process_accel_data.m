function process_accel_data(filename, varargin)

opt.burstparams = 'ask';  % or 'file','gui'
opt = parsevarargin(opt,varargin,2);

accelnames = {'FwdAcc','SideAcc','UpAcc'};
orientnames = {'Roll','Pitch','Yaw'};
kinnames = {'FwdPos','FwdVel','FwdAccVid','SidePos','SideVel','SideAccVid','TailFwd','TailSide'};

if (nargin == 0) || isempty(filename)
    [fn,pn] = uigetfile('*.mat','Choose full data file');
    filename = fullfile(pn,fn);
end

F = load(filename);
varnames = fieldnames(F);

emgnames = varnames(~ismember(varnames,[accelnames orientnames kinnames 't']));

t = F.t;
N = length(t);
emg = zeros(N,length(emgnames));
for i = 1:length(emgnames)
    emg(:,i) = F.(emgnames{i});
end
nemg = size(emg,2);

%first find peaks in tail motion
[pk,pkind,pksgn] = findpeaks2(F.TailSide,'minmax','sort','none');
if (pksgn(1) == -1)
    pk = pk(2:end);
    pkind = pkind(2:end);
    pksgn = pksgn(2:end);
end
if (pksgn(end) == 1)
    pkind = pkind(1:end-1);
    pk = pk(1:end-1);
    pksgn = pksgn(1:end-1);
end
pkt = t(pkind);

phase1 = (0:length(pk)-1)/2;
tailphase = zeros(size(t));
tailphase(pkind(1):pkind(end)) = interp1(pkt,phase1, t(pkind(1):pkind(end)), 'spline');

if strcmp(opt.burstparams,'ask')
    fprintf('Get burst identification parameters from\n  (F) file\n  (G) GUI\n');
    c = input('? ','s');
    switch lower(c)
        case 'f'
            opt.burstparams = 'file';
        case 'g'
            opt.burstparams = 'gui';
        otherwise
            opt.burstparams = 'gui';
    end
end

switch lower(opt.burstparams)
    case 'gui'
        emgdata = findbursts_gui(t,emg);
        
    otherwise
        if strcmp(opt.burstparams,'file')
            [pn,fn,ext] = fileparts(filename);
            [fn,pn] = uigetfile(fullfile(pn,'*.mat'),'Choose data file with threshold');
            burstidfile = fullfile(pn,fn);
        elseif exist(opt.burstparams,'file')
            burstidfile = opt.burstparams;
        end
        B = load(burstidfile,'interburstdur','spikethreshold','minspikes');
        if ~isfield(B,'interburstdur') || ~isfield(B,'spikethreshold') || ~isfield(B,'minspikes')
            error('No burst info in data file');
        end
        
        emgdata = findbursts_gui(t,emg,'interburstdur',B.interburstdur, ...
            'minspikes',B.minspikes,'threshold',B.spikethreshold, 'quiet');
end

spiket = emgdata.spiket;
spikeamp = emgdata.spikeamp;
spikethreshold = emgdata.spikethreshold;
interburstdur = [emgdata.burst.interburstdur];
minspikes = [emgdata.burst.minspikes];

burston0 = emgdata.burston;
burstoff0 = emgdata.burstoff;

burstonphase0 = interp1(t,tailphase, burston0);
burstoffphase0 = interp1(t,tailphase, burstoff0);

burstint0 = zeros(size(burston0));
for i = 1:size(burston0,1)
    for c = 1:nemg,
        isburst = (t >= burston0(i,c)) & (t <= burstoff0(i,c));
        burstint0(i,c) = trapz(t(isburst),abs(emg(isburst,c)));
        burstamp0(i,c) = nanmean(abs(emg(isburst,c)));
    end
end

nbeats = length(pkind)/2;
tbeat = zeros(nbeats-1,1);
acc = zeros(nbeats-1,1);
accpk = zeros(nbeats-1,1);
per = zeros(nbeats-1,1);
vel = zeros(nbeats-1,1);
tailamp = zeros(nbeats-1,1);
burstint = NaN(nbeats-1,nemg);
burstamp = NaN(nbeats-1,nemg);
burstonphase = NaN(nbeats-1,nemg);
burstoffphase = NaN(nbeats-1,nemg);
burstduty = NaN(nbeats-1,nemg);
nbursts = zeros(nbeats-1,nemg);
for i = 1:nbeats-1
    j = 2*i - 1;
    isbeat = pkind(j):pkind(j+2);
    
    tbeat(i) = pkt(j);
    
    acc(i) = nanmean(F.FwdAcc(isbeat)) * 1000;
    accpk(i) = prctile(F.FwdAcc(isbeat),99) * 1000;
    per(i) = pkt(j+2) - pkt(j);
    vel(i) = nanmean(F.FwdVel(isbeat));
    tailamp(i) = ((pk(j) - pk(j+1)) + (pk(j+2) - pk(j+1)))/2;
    
    for k = 1:nemg
        burstind = find((burston0(:,k) >= pkt(j)) & (burston0(:,k) < pkt(j+2)));
        nbursts(i,k) = length(burstind);
        if (length(burstind) > 1)
            dur1 = (burstoff0(burstind,k) - burston0(burstind,k));
            [~,longest] = max(dur1);
            burstind = burstind(longest);
        end
        
        if ~isempty(burstind)
            ph0 = tailphase(pkind(j));
            burstonphase(i,k) = burstonphase0(burstind,k) - ph0;
            burstoffphase(i,k) = burstoffphase0(burstind,k) - ph0;
            burstint(i,k) = burstint0(burstind,k);
            burstamp(i,k) = burstamp0(burstind,k);
            
            burstduty(i,k) = (burstoff0(burstind,k) - burston0(burstind,k)) / per(i);
        end
    end        
end

[pn,fn,ext] = fileparts(filename);
outname = fullfile(pn,[fn '-analysis.mat']);

[fn,pn] = uiputfile('*.mat','Choose output file',outname);
outname = fullfile(pn,fn);

save(outname,'spiket','spikeamp','spikethreshold','interburstdur','minspikes',...
    'burston0', 'burstoff0', 'tbeat','acc','accpk','vel','per','burstint',...
    'burstonphase','burstoffphase', 'burstduty','burstamp','nbursts','tailamp','emgnames');

colnames = {'tbeat','per','vel','accmn','accpk','tailamp'};
burstnames = {'int','amp',...
    'onph','offph','duty','nincycle'};
burstnames = repmat(burstnames,[length(emgnames) 1]);
for i = 1:length(emgnames)
    for j = 1:length(burstnames)
        burstnames{i,j} = [burstnames{i,j} '-' emgnames{i}];
    end
end
colnames = [colnames burstnames(:)'];
colunits = {'s','s','mm/s','mm/s^2','mm/s^2','mm'};
burstunits = repmat({'V*s','V','','','',''},[length(emgnames) 1]);
colunits = [colunits burstunits(:)'];

tplt1 = repmat('%s,',[1 length(colnames)]);
tplt1 = [tplt1(1:end-1) '\n'];

[pn,fn,ext] = fileparts(outname);
outname2 = fullfile(pn,[fn '.csv']);

fid = fopen(outname2,'w');
fprintf(fid,tplt1,colnames{:});
fprintf(fid,tplt1,colunits{:});

tplt2 = repmat('%.6f,',[1 length(colnames)]);
tplt2 = [tplt2(1:end-1) '\n'];

X = [tbeat,per,vel,acc,accpk,tailamp,burstint,burstamp,burstonphase,burstoffphase,burstduty,nbursts];
assert(size(X,2) == length(colnames));

fprintf(fid,tplt2,X');
fclose(fid);

    






