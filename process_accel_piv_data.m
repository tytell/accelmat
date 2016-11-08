function P = process_accel_piv_data(filename, kin, varargin)

opt.nframesmean = 5;
opt.leftsidevortexsign = -1;
opt.showdiagnostics = true;
opt.sampfreq = [];      % if empty, use the same as the kinematics data
opt.frameskip = 1;
opt.tend = 0;           % last time in the video. Usually <0 if a post trigger was used

opt = parsevarargin(opt,varargin, 3);

if isempty(opt.sampfreq)
    dt = kin.t(2) - kin.t(1);
else
    dt = 1/opt.sampfreq;
end

F = load(filename);

names = fieldnames(F);

tok = regexp(names,'Cl(\d+)circ((N|ctr)?)','tokens','once');

ismatch = ~cellfun(@isempty,tok);
num = zeros(size(ismatch));
num(ismatch) = cellfun(@(x) str2double(x{1}), tok(ismatch));
var = cell(size(ismatch));
var(ismatch) = cellfun(@(x) x{2}, tok(ismatch), 'UniformOutput',false);
[var{~ismatch}] = deal('');

iscirc = false(size(names));
iscirc(ismatch) = cellfun(@(x) isempty(x{2}), tok(ismatch));

indcirc = find(iscirc);
indcirc = indcirc(:)';
if isempty(indcirc)
    error('Tyler - this is a weird file.  Do something!');
end

nframes = length(F.(names{indcirc(1)}));
dur = nframes * dt;

frameskip = opt.frameskip;
fr = (0:nframes-1)*frameskip + 1;
t = fr*dt;
t = t - t(end) + opt.tend;

if range(t) > 1.05*range(kin.t)
    warning('PIV data is longer than kinematic data. Only analyzing overlapping times.');
    nframes = length(kin.t);
elseif range(kin.t) > 1.05*range(t)
%     fprintf('Kinematics time range = [%f %f]; PIV time range = [%f %f]\n', ...
%         kin.t(1), kin.t(end), 
    warning('Kinematics data is longer than PIV data.  Only analyzing overlapping times.');
    frameskip = round(length(kin.t) / nframes);
    
    fprintf('PIV Frameskip seems to be %d.\n', frameskip);
    fs1 = input(sprintf('What is the correct frameskip (default = %d)? ', frameskip));
    if ~isempty(fs1) && (fs1 >= 1)
        frameskip = fs1;
    end
end

span = (t >= kin.t(1)) & (t <= kin.t(end));
trange = [first(t, span), last(t, span)];
nframes = sum(span);
fr = fr(span);

circ0 = NaN(nframes,length(indcirc));

for i = 1:length(indcirc)
    j = indcirc(i);
    good = ~cellfun(@isempty, F.(names{j})) & span;
    
    %make sure we only take overlapping times
    good2 = good(span);
    
    circ0(good2,num(j)) = cat(1,F.(names{j}){good});
end
circ0(circ0 == 0) = NaN;

indctr = find(ismember(var,{'ctr'}));
vxx0 = NaN(nframes,length(indcirc));
vxy0 = NaN(nframes,length(indcirc));

for i = 1:length(indctr)
    j = indctr(i);
    good = ~cellfun(@isempty, F.(names{j})) & span;
    
    %make sure we only take up through nframes frames
    good2 = good(span);

    if all(~good)
        continue;
    end
    
    xy = cat(1,F.(names{j}){good});
    
    vxx0(good2,num(j)) = xy(:,1);
    vxy0(good2,num(j)) = xy(:,2);
end
bad = (vxx0 == 0) & (vxy0 == 0);
vxx0(bad) = NaN;
vxy0(bad) = NaN;

%get rid of any elements with no circulation values
good = any(isfinite(circ0));
circ0 = circ0(:,good);
vxx0 = vxx0(:,good);
vxy0 = vxy0(:,good);
nvx = size(circ0,2);

%now match sampling frequencies
if opt.sampfreq / frameskip > 1/(kin.t(2) - kin.t(1))
    %upsample the kinematics data
    tpeak = NaN(size(kin.indpeak));
    good = isfinite(kin.indpeak);
    tpeak(good) = kin.t(kin.indpeak(good));
    
    indpeak = round((tpeak - t(1)) * opt.sampfreq)+1;
    
    swimvecx = NaN(size(t));
    swimvecy = NaN(size(t));
    good = isfinite(kin.swimvecx);
    swimvecx(span) = interp1(kin.t(good),kin.swimvecx(good), t(span));
    swimvecy(span) = interp1(kin.t(good),kin.swimvecy(good), t(span));   
else
    indpeak = kin.indpeak;
    swimvecx = kin.swimvecx;
    swimvecy = kin.swimvecy;
end
    
circ1 = NaN(max(fr), nvx);
circ1(fr,:) = circ0;
vxx1 = NaN(max(fr), nvx);
vxx1(fr,:) = vxx0;
vxy1 = NaN(max(fr), nvx);
vxy1(fr,:) = vxy0;

circ0 = circ1;
vxx0 = vxx1;
vxy0 = vxy1;

vxt0 = repmat(t(:),[1 size(circ0,2)]);
vxt0(isnan(circ0)) = NaN;
if (size(vxt0,1) > size(circ0,1))
    vxt0 = vxt0(1:size(circ0,1),:);
end

vxstart = first(isfinite(circ0));
vxsgn = nanmedian(sign(circ0));

ind = indpeak(end,:);

tailside = kin.sidepeak(end,:);
ind2 = cat(2,ind,nframes);

iscircsignmatch = false(length(ind),nvx);
circmean = zeros(length(ind),nvx);
fracdef = zeros(length(ind),nvx);
for i = 1:length(ind)
    if (ind2(i) < 1) || (ind2(i+1) < 1) || ...
            (ind2(i) > length(t)) || (ind2(i+1) > length(t))
        warning('Kinematics tail beat %i is out of the PIV time range. Skipping', i);
        continue;
    end
    
    if tailside(i) == 'L'
        findsign = opt.leftsidevortexsign;
    else
        findsign = -opt.leftsidevortexsign;
    end  
%     %look for matching sign of the circulation
%     circsignmatch1 = sign(circ0(ind(i):ind(i+1),:)) == findsign;
%     %call it a match if 90% of frames have the right circulation
%     iscircsignmatch(i,:) = sum(circsignmatch1) > (0.9 * (ind(i+1)-ind(i)+1));
    
    circmean(i,:) = nanmean(circ0(ind2(i):ind2(i+1),:));
    iscircsignmatch(i,:) = sign(circmean(i,:)) == findsign;
    fracdef(i,:) = sum(isfinite(circ0(ind2(i):ind2(i+1),:))) ./ (ind2(i+1)-ind2(i)+1);
end

%match vortices to tailbeats by looking for the vortex that's present with
%the right sign in the largest fraction of frames during each tailbeat,
%making sure that we don't pick the same vortex twice.  Assumes that
%vortices are matched in order (so that if vortex 5 is matched to tailbeat
%4, vortex 6 couldn't match tailbeat 3)

circmatch = fracdef .* double(iscircsignmatch);
vxmatch = NaN(length(ind),1);
a = 1;
for i = 1:size(circmatch,1)
    %first look for ones that are at least 80% in this cycle
    k = find(circmatch(i,a:end) >= 0.8/frameskip);
    k = k + a-1;
    if length(k) == 1
        vxmatch1 = k;
    elseif length(k) > 1
        [~,vxmatch1] = max(abs(circmean(i,k)));
        vxmatch1 = k(vxmatch1);
    else
        [v1,vxmatch1] = max(circmatch(i,a:end));
        vxmatch1 = vxmatch1+a-1;
        if v1 <= 0.5
            vxmatch1 = [];
        end
    end
    
    if ~isempty(vxmatch1)
        vxmatch(i) = vxmatch1;
        a = vxmatch1+1;
    end
end

vxcirc = NaN(size(ind));
vxt = NaN(size(ind));
vxprevcirc = NaN(size(ind));
vxdist = NaN(size(ind));
vxcircstd = NaN(size(ind));
vxprevcircstd = NaN(size(ind));
vxdiststd = NaN(size(ind));
vxx = NaN(size(ind));
vxy = NaN(size(ind));
vxprevx = NaN(size(ind));
vxprevy = NaN(size(ind));
vxang = NaN(size(ind));
vxind = NaN(size(ind));

for i = 1:length(ind)
    if ~isfinite(vxmatch(i))
        warning('No vortex found for tail beat %d',i);
        continue
    end
    
    vx1 = vxmatch(i);
    
    %start at the beginning of the vortex
    k1 = vxstart(vx1);
    if (k1 < ind(i))
        %unless the beginning of the vortex is before the beginning of this
        %tail beat
        k1 = ind(i);
    end
    
    k = k1+(0:(opt.nframesmean-1)*frameskip);
    k = k((k >= 1) & (k <= size(circ0,1)));
    
    vxcirc(i) = nanmean(circ0(k,vx1));
    vxt(i) = nanmean(vxt0(k,vx1));
    vxcircstd(i) = nanstd(circ0(k,vx1));
    vxx1 = vxx0(k,vx1);
    vxy1 = vxy0(k,vx1);
    vxx(i) = vxx1(1);
    vxy(i) = vxy1(1);
    
    vxprev1 = last((sum(isfinite(circ0(k,1:vx1-1))) >= (length(k)-2)/frameskip) & ...
        ~iscircsignmatch(i,1:vx1-1));
    if ~isempty(vxprev1)
        k = vxstart(vxprev1)+(0:opt.nframesmean-1);
        k = k((k >= 1) & (k <= size(circ0,1)));

        vxprevcirc(i) = nanmean(circ0(k,vxprev1));
        vxprevcircstd(i) = nanstd(circ0(k,vxprev1));

        if length(k) > length(vxx1)
            k = k(1:length(vxx1));
        end
        vxxprev1 = vxx0(k,vxprev1);
        vxyprev1 = vxy0(k,vxprev1);
        vxprevx(i) = vxxprev1(1);
        vxprevy(i) = vxyprev1(1);
        
        vxdist1 = sqrt((vxx1 - vxxprev1).^2 + (vxy1 - vxyprev1).^2);
        
        vxdist(i) = nanmean(vxdist1);
        vxdiststd(i) = nanstd(vxdist1);
        
        %previous vortex should always have > x than the current one,
        %because it's further down in the wake
        dx = vxxprev1 - vxx1;
        dy = vxyprev1 - vxy1;
        
        %angle of the vortex pair to the opposite of the swimming direction
        %(ie, backwards)
        dax = -dx.*swimvecx(k)' - dy.*swimvecy(k)';
        dlat = +dx.*swimvecy(k)' - dy.*swimvecx(k)';
        
        vxang1 = atan2(dlat,dax);
        vxang(i) = angmean(vxang1);
    end
end

if opt.showdiagnostics
    tailpos = kin.exc(end,:);
    tailpos = (tailpos - nanmean(tailpos)) ./ range(tailpos);
    
    good = isfinite(kin.indpeak(end,:));
    tailpeak = tailpos(kin.indpeak(end,good));
    tpeak = kin.tpeak(end,good);
    
    circs = circ0 ./ range(circ0(:));
    vxcircs = vxcirc ./ range(circ0(:));
    
    plot(kin.t, tailpos,'r-', tpeak,tailpeak,'ro');
    
    goodcirc = any(isfinite(circs),2);
    addplot(vxt0(goodcirc,:), circs(goodcirc,:), 'k-');
    
    addplot(cat(1,tpeak,vxt), cat(1,tailpeak,vxcircs), 'b*-');
    
    notmatched = isnan(vxt);
    addplot(tpeak(notmatched), tailpeak(notmatched), 'gs');
    drawnow;
    if ~inputyn('Continue?','default',true)
        error('Aborted...');
    end
end

P.vxx = vxx;
P.vxy = vxy;
P.vxprevx = vxprevx;
P.vxprevy = vxprevy;
P.vxang = vxang;
P.vxcirc = vxcirc;
P.vxcircstd = vxcircstd;
P.vxprevcirc = vxprevcirc;
P.vxprevcircstd = vxprevcircstd;
P.vxdist = vxdist;
P.vxdiststd = vxdiststd;
