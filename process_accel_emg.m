function emg = process_accel_emg(emg, kin, varargin)

opt.interburstdur = [];
opt.spikethreshold = [];
opt.minspikes = [];
opt.goodchan = logical([]);
opt.override = struct([]);
opt.trigger = 'falling';
opt.emgposition = struct([]);
opt.doplot = false;

opt = parsevarargin(opt, varargin, 3);

t = emg.t;
chan = struct2cell(emg.channels);
chan = cat(2,chan{:});
channames = fieldnames(emg.channels);
nchan = length(channames);

trigind = find(strcmpi(channames,'trigger'));
if length(trigind) == 1
    trig = chan(:,trigind);
    switch opt.trigger
        case 'falling'
            t0 = first(t(1:end-1),(trig(1:end-1) > 3) & (trig(2:end) <= 3));
        case 'rising'
            t0 = first(t(1:end-1),(trig(1:end-1) < 3) & (trig(2:end) >= 3));
    end            
    chan = chan(:,[1:trigind-1 trigind+1:end]);
    channames = channames([1:trigind-1 trigind+1:end]);
    nchan = nchan-1;
else
    warning('Could not find trigger channel.  Assuming trigger at end of data');
    t0 = t(end);
end
t = t - t0;

if isempty(opt.emgposition)
    error('Need to have positions for the EMG electrodes');
end

emgposition = zeros(nchan,2);
for i = 1:length(channames)
    if ~isfield(opt.emgposition, channames{i})
        error('Could not find a position for channel %s in emgpositions.',channames{i});
    end
    emgposition(i,:) = opt.emgposition.(channames{i});
end

span = find((t >= kin.t(1)) & (t <= kin.t(end)));
if ~isempty(opt.spikethreshold) && ~isempty(opt.interburstdur) && ...
        ~isempty(opt.minspikes) && ~isempty(opt.goodchan)
    data = findbursts_gui(t(span), chan(span,:), 'interburstdur',opt.interburstdur, ...
        'threshold',opt.spikethreshold, 'minspikes',opt.minspikes, ...
        'goodchan',opt.goodchan, 'override',opt.override);
else
    data = findbursts_gui(t(span),chan(span,:));
end

nburst = size(data.burston,1);
nwave = size(kin.tcurvepeak,2);

%interpolate curvature at each emg position
nfr = size(kin.curve,2);
emgcurve = NaN(nchan,nfr);
for f = 1:size(kin.curve,2)
    s1 = kin.s;
    c1 = kin.curve(:,f);
    c1([1 end]) = 0;
    emgcurve1 = interp1(s1,c1, emgposition(:,2));
    emgcurve(:,f) = emgcurve1;    
end

burstdur0 = (data.burstoff - data.burston)';
burstfreq0 = NaN(size(data.burstt));
for c = 1:nchan
    good = isfinite(data.burstt(:,c));
    burstt1 = data.burstt(good,c);
    freq1 = NaN(size(burstt1));
    freq1(2:end-1) = 2./(burstt1(3:end) - burstt1(1:end-2));
    
    burstfreq0(good,c) = freq1;
end
burstfreq0 = burstfreq0';

burstamp0 = NaN(nchan,nburst);
burstint0 = NaN(nchan,nburst);
for b = 1:nburst
    for c = 1:nchan
        if (isfinite(data.burston(b,c)) && isfinite(data.burstoff(b,c)))
            isburst = (t >= data.burston(b,c)) & (t <= data.burstoff(b,c));
            burstamp0(c,b) = nanmean(abs(chan(isburst,c)));
            burstint0(c,b) = trapz(t(isburst),abs(chan(isburst,c)));
        end
    end
end

%and find the phase
emgcurvephase = NaN(nchan,nfr);
burstonphase = NaN(nchan,nwave);
burstoffphase = NaN(nchan,nwave);
burston = NaN(nchan,nwave);
burstoff = NaN(nchan,nwave);
burstctr = NaN(nchan,nwave);
burstdur = NaN(nchan,nwave);
burstamp = NaN(nchan,nwave);
burstint = NaN(nchan,nwave);
burstfreq = NaN(nchan,nwave);

for c = 1:nchan
    if emgposition(c,1) == 2
        sgn = -1;
    else
        sgn = 1;
    end
    
    [pk,pkind,pksgn] = findpeaks2(sgn*emgcurve(c,:), 'minmax');
    pkt = kin.t(pkind);
    
    %construct a phase based on the peaks in curvature
    pkphase1 = zeros(size(pk));
    pkphase1(pksgn == -1) = 0.5;
    
    pkcycle1 = zeros(size(pkt));
    pkcycle1(pksgn == 1) = 1;
    pkcycle1 = cumsum(pkcycle1);
    pkphase1 = pkphase1 + pkcycle1;
    
    if (min(pkphase1) < 1)
        pkphase1 = pkphase1 - floor(min(pkphase1)) + 1;
    end
    
    iscurvet = (kin.t >= min(pkt)) & (kin.t <= max(pkt));
    
    curvephase1 = NaN(size(kin.t));
    curvephase1(iscurvet) = interp1(pkt,pkphase1, kin.t(iscurvet), 'spline');
    
    pkt = pkt(pksgn == 1);
    
    %figure out which identified wave the new curvature peak corresponds to
    iscurvesgn = mode(sign(kin.ampcurvepeak)) == sgn;
    wavet = min(kin.tcurvepeak);
    wavet(~iscurvesgn) = NaN;
    pkwavenum = NaN(size(pkt));
    for i = 1:length(pkt)
        d = pkt(i)-wavet;
        d(d < 0) = NaN;
        [~,ind] = min(d);
        pkwavenum(i) = ind;
    end
       
    emgcurvephase(c,:) = curvephase1;
    
    %and put the burst phases in a matrix of the same size
    burstonphase1 = interp1(kin.t,curvephase1, data.burston(:,c));
    burstoncycle1 = floor(burstonphase1);
    burstoffphase1 = interp1(kin.t,curvephase1, data.burstoff(:,c));
    
    burstonphase2 = NaN(1,nwave);
    burstoffphase2 = NaN(1,nwave);
    good = isfinite(burstoncycle1) & (burstoncycle1 >= 1) & (burstoncycle1 <= length(pkwavenum));
    
    burstwavenum = pkwavenum(burstoncycle1(good));
    burstonphase2(burstwavenum) = burstonphase1(good) - burstoncycle1(good);
    burstoffphase2(burstwavenum) = burstoffphase1(good) - burstoncycle1(good);
    
    burstonphase(c,:) = burstonphase2;
    burstoffphase(c,:) = burstoffphase2;
    
    %get the actual burst times
    burston(c,burstwavenum) = data.burston(good,c)';
    burstoff(c,burstwavenum) = data.burstoff(good,c)';
    burstctr(c,burstwavenum) = data.burstt(good,c)';

    burstfreq(c,burstwavenum) = burstfreq0(c,good);
    
    %also get the duty cycle and burst amplitude and intensity
    burstdur(c,burstwavenum) = burstdur0(c,good);
    burstamp(c,burstwavenum) = burstamp0(c,good);
    burstint(c,burstwavenum) = burstint0(c,good);
end

burstduty = burstdur ./ repmat(kin.per,[nchan 1]);

if opt.doplot
    emgpos = repmat(emgposition(:,2),[1 nwave]);
    emgside = repmat(emgposition(:,1),[1 nwave]);

    jit1 = rand(size(burstonphase)) * 0.02;
    
    isleft = emgside == 1;
    X = emgpos(isleft)+jit1(isleft)-0.011;
    plot([X(:) X(:)]', [burstonphase(isleft) burstoffphase(isleft)]', 'bo-');
    
    X = emgpos(~isleft)+jit1(~isleft)+0.011;
    addplot([X(:) X(:)]', [burstonphase(~isleft) burstoffphase(~isleft)]', 'gs-');
    
    [~,~,emgind] = unique(emgpos);
    onmn = accumarray(emgind(:),burstonphase(:), [], @(x) angmean(2*pi*x));
    onmn = onmn/(2*pi);
    offmn = accumarray(emgind(:),burstoffphase(:), [], @(x) angmean(2*pi*x));
    offmn = offmn/(2*pi);
end

emg.pos = emgposition(:,2);
emg.side = emgposition(:,1);
emg.curve = emgcurve;
emg.burst = data.burst;
emg.spikethreshold = data.spikethreshold;
emg.interburstdur = data.interburstdur;
emg.minspikes = data.minspikes;
emg.goodchan = data.goodchan;
emg.burstoverride = data.override;
emg.burston = burston;
emg.burstoff = burstoff;
emg.burstctr = burstctr;
emg.burstfreq = burstfreq;
emg.curvephase = emgcurvephase;
emg.burstonphase = burstonphase;
emg.burstoffphase = burstoffphase;
emg.burstduty = burstduty;
emg.burstdur = burstdur;
emg.burstamp = burstamp;
emg.burstint = burstint;



