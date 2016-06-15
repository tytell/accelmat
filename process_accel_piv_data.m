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
if isempty(indcirc)
    error('Tyler - this is a weird file.  Do something!');
end

nframes = length(F.(names{indcirc(1)}));

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

ind(end+1:end+2) = size(circ0,1);
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
    
    vxprev1 = last((vxstart < ind(i)) & (vxsgn == -vxsgn(vx1)));
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
        
        vxdist1 = sqrt((vxx1 - vxxprev1).^2 + (vxy1 - vxyprev1).^2);
        
        vxdist(i) = nanmean(vxdist1);
        vxdiststd(i) = nanstd(vxdist1);
    end
end

P.vxcirc = vxcirc;
P.vxcircstd = vxcircstd;
P.vxprevcirc = vxprevcirc;
P.vxprevcircstd = vxprevcircstd;
P.vxdist = vxdist;
P.vxdiststd = vxdiststd;
