function P = process_accel_piv_data(filename, kin, varargin)

opt.nframesmean = 5;
opt.leftsidevortexsign = -1;

opt = parsevarargin(opt,varargin, 3);

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

nframes = length(F.(names{indcirc(1)}));
if nframes > length(kin.t)
    warning('PIV data has more frames than kinematic data. Only analyzing overlapping frames.');
    nframes = length(kin.t);
elseif length(kin.t) > nframes+5
    warning('Kinematics data has more frames than PIV data.  Only analyzing overlapping frames.');
end

circ0 = NaN(nframes,length(indcirc));

for i = 1:length(indcirc)
    j = indcirc(i);
    good = ~cellfun(@isempty, F.(names{j}));
    circ0(good,num(j)) = cat(1,F.(names{j}){good});
end

indctr = find(ismember(var,{'ctr'}));
vxx0 = NaN(nframes,length(indcirc));
vxy0 = NaN(nframes,length(indcirc));

for i = 1:length(indctr)
    j = indctr(i);
    good = ~cellfun(@isempty, F.(names{j}));
    
    xy = cat(1,F.(names{j}){good});
    
    vxx0(good,num(j)) = xy(:,1);
    vxy0(good,num(j)) = xy(:,2);
end

vxstart = first(isfinite(circ0));
vxsgn = nanmedian(sign(circ0));
circ0mn = nanmean(circ0);

ind = kin.indpeak(end,:);
tailside = kin.sidepeak(end,:);

vxcirc = NaN(size(ind));
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

ind(end+1:end+2) = nframes;
for i = 1:length(ind)-2
    if tailside(i) == 'L'
        findsign = opt.leftsidevortexsign;
    else
        findsign = -opt.leftsidevortexsign;
    end  
    
    indvortex = find((vxstart >= ind(i)) & (vxstart < ind(i+2)) & ...
        (vxsgn == findsign));
    [~,vx1] = max(circ0mn(indvortex));
    vx1 = indvortex(vx1);

    if isempty(vx1) || ((i < length(ind)-1) && (vx1 >= ind(i+2)))
        warning('No vortex found for tail beat %d\n',i);
        continue;
    end
    
    k = vxstart(vx1)+(0:opt.nframesmean-1);
    k = k((k >= 1) & (k <= size(circ0,1)));
    vxcirc(i) = nanmean(circ0(k,vx1));
    vxcircstd(i) = nanstd(circ0(k,vx1));
    vxx1 = vxx0(k,vx1);
    vxy1 = vxy0(k,vx1);
    vxx(i) = vxx1(1);
    vxy(i) = vxy1(1);
    
    vxprev1 = last(all(isfinite(circ0(k,1:vx1-1))) & (vxsgn(1:vx1-1) == -vxsgn(vx1)));
    if ~isempty(vxprev1)
        vxprevcirc(i) = nanmean(circ0(k,vxprev1));
        vxprevcircstd(i) = nanstd(circ0(k,vxprev1));

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
        dax = -dx.*kin.swimvecx(k)' - dy.*kin.swimvecy(k)';
        dlat = +dx.*kin.swimvecy(k)' - dy.*kin.swimvecx(k)';
        
        vxang1 = atan2(dlat,dax);
        vxang(i) = angmean(vxang1);
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
