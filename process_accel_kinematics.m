function K = process_accel_kinematics(filename, varargin)

opt.view = 'ventral';
opt.mirror = false;     % filmed through a mirror at 45deg (flips y axis)

opt = parsevarargin(opt,varargin, 2);

load(filename,'exs','eys','scale','fishlenmm','haxmmss','haymmss','humms','hvmms',...
    'hxmm','hymm','hxs','hys','t','txmm','tymm','txs','tys');

%reorder the points so that they go from head to tail
%first get the head to tail axis
htx = txs - hxs;
hty = tys - hys;
htmag = sqrt(htx.^2 + hty.^2);
htx = htx./htmag;
hty = hty./htmag;

nextra = size(exs,1);

%then project the extra points on to the head to tail axis
htdist = (exs - repmat(hxs,[nextra 1])) .* repmat(htx,[nextra 1]) + ...
    (eys - repmat(hys,[nextra 1])) .* repmat(hty,[nextra 1]);

[~,ord] = sort(htdist);
good = all(isfinite(htdist));
k = first(good);
if all(all(ord(:,good) == repmat(ord(:,k),[1 sum(good)])))
    ord = ord(:,k);
    exs = exs(ord,:);
    eys = eys(ord,:);
else
    for i = 1:size(ord,2)
        exs(:,i) = exs(ord(:,i),i);
        eys(:,i) = eys(ord(:,i),i);
    end
end

mx = [hxs; exs; txs];
my = [hys; eys; tys];
npts = size(mx,1);

if opt.mirror
    my = -my;
end

mxmm = mx*scale;
mymm = my*scale;
dsmm = sqrt(diff(mxmm).^2 + diff(mymm).^2);
smm = [zeros(1,size(mx,2)); cumsum(dsmm)];

swimvecx = NaN(1,size(mxmm,2));
swimvecy = NaN(1,size(mxmm,2));
midx = NaN(size(mxmm));
midy = NaN(size(mymm));

ind = find(all(isfinite(mxmm)));
for i = ind
    %straight line interpolant, minimizing the perpendicular distance
    %to the points (doesn't prioritize x or y, and works if we're
    %horizontal or vertical)
    sp1 = spap2(1, 2, smm(:,i)', [mxmm(:,i),mymm(:,i)]');
    xy1 = fnval(sp1,smm(:,i)')';
    midx(:,i) = xy1(:,1);
    midy(:,i) = xy1(:,2);
    
    swimvecx1 = xy1(end,1) - xy1(1,1);
    swimvecy1 = xy1(end,2) - xy1(1,2);
    mag = sqrt(swimvecx1.^2 + swimvecy1.^2);
    
    swimvecx(i) = swimvecx1/mag;
    swimvecy(i) = swimvecy1/mag;
end

excx = mxmm - midx;
excy = mymm - midy;

switch opt.view
    case 'ventral'
        leftsign = -1;
    case 'dorsal'
        leftsign = 1;
end

exc = -excx.*repmat(swimvecy,[size(excx,1) 1]) + excy.*repmat(swimvecx,[size(excx,1) 1]);

pkind = cell(1,2);
pkamp = cell(1,2);
pksite = cell(1,2);
for p = 1:npts
    [pkamp1,pkind1,pksgn1] = findpeaks2(exc(p,:),'minmax');
    
    isleft = pksgn1 == leftsign;
    pkind{1} = cat(2,pkind{1},pkind1(isleft));
    pkamp{1} = cat(2,pkamp{1},pkamp1(isleft));
    pksite{1} = cat(2,pksite{1},p * ones(1,sum(isleft)));
    
    pkind{2} = cat(2,pkind{2},pkind1(~isleft));
    pkamp{2} = cat(2,pkamp{2},pkamp1(~isleft));
    pksite{2} = cat(2,pksite{2},p * ones(1,sum(~isleft)));
end

for side = 1:2
    [pkind{side},ord] = sort(pkind{side});
    pkamp{side} = pkamp{side}(ord);
    pksite{side} = pksite{side}(ord);
end

npk = sum(pksite{1} == 1) + sum(pksite{1} == 1) + 3;

indpeak = NaN(npts,npk);
amppeak = NaN(npts,npk);
sidepeak = spaces(npts,npk);
for side = 1:2
    pkn = side;
    prevsite = pksite{side}(1)-1;
    if side == 1
        sidestr = 'L';
    else
        sidestr = 'R';
    end
    pkind1 = pkind{side};
    pkamp1 = pkamp{side};
    pksite1 = pksite{side};
    
    for a = 1:length(pkind1)
        site = pksite1(a);
        if (site == prevsite + 1)
            if isnan(indpeak(site,pkn)) || (amppeak(site,pkn) < pkamp1(a))
                indpeak(site,pkn) = pkind1(a);
                amppeak(site,pkn) = pkamp1(a);
                sidepeak(site,pkn) = sidestr;
            end
            prevsite = site;
        elseif (site == prevsite - 1) && ...
                (abs(pkind1(a) - indpeak(prevsite,pkn)) < 10)
            if isnan(indpeak(site,pkn)) || (amppeak(site,pkn) < pkamp1(a))
                indpeak(site,pkn) = pkind1(a);
                amppeak(site,pkn) = pkamp1(a);
                sidepeak(site,pkn) = sidestr;
            end
            %don't update prevsite, because we skipped back in the order to
            %match this cycle
        elseif (site < prevsite)
            pkn = pkn+2;
            indpeak(site,pkn) = pkind1(a);
            amppeak(site,pkn) = pkamp1(a);
            sidepeak(site,pkn) = sidestr;
            prevsite = site;
        elseif (pkn > side) && (site > 2) && isfinite(indpeak(site-1,pkn-2)) && ...
                isfinite(indpeak(site-2,pkn-2))
            pkdiff1 = indpeak(site-1,pkn-2)-indpeak(site-2,pkn-2);
            pkdiff2 = pkind1(a)-indpeak(site-1,pkn-2);
            
            if (pkdiff2 < pkdiff1/dsmm(site-2,indpeak(site-2,pkn-2))*fishlenmm)
                indpeak(site,pkn-2) = pkind1(a);
                amppeak(site,pkn-2) = pkamp1(a);
                sidepeak(site,pkn-2) = sidestr;
            else
                warning('Could not match a peak!');
            end
            %don't update prevsite, because we had to skip back a cycle to
            %match this peak
        end
    end
end

good = any(isfinite(indpeak));
indpeak = indpeak(:,good);
amppeak = amppeak(:,good);
sidepeak = sidepeak(:,good);

[~,ord] = sort(nanmean(indpeak));
indpeak = indpeak(:,ord);
amppeak = amppeak(:,ord);
sidepeak = sidepeak(:,ord);

ang = atan2(diff(mymm),diff(mxmm));
dang = diff(ang);
curve = NaN(size(mymm));
curve(2:end-1,:) = dang./((dsmm(1:end-1,:) + dsmm(2:end,:))/2);

%positive curvature should be concave toward the right
curve = curve * (-leftsign);

%cheating a bit - subtract off the mean curvature at each point so that
%curvature is centered around zero
curve = curve - repmat(nanmean(curve,2),[1 size(curve,2)]);

indcurvepeak = NaN(size(indpeak));
ampcurvepeak = NaN(size(indpeak));
for pt = 2:npts-1
    for pk = 1:size(indpeak,2)-1
        if isfinite(indpeak(pt,pk)) && isfinite(indpeak(pt,pk+1)) && ...
                indpeak(pt,pk+1) > indpeak(pt,pk)
            pkind1 = indpeak(pt,pk):indpeak(pt,pk+1)-1;
            curve1 = curve(pt,pkind1);
            if sidepeak(pt,pk) == 'L'
                sidemult = -1;
            else
                sidemult = 1;
            end
            [maxcurve1,indcurve1] = max(sidemult*curve1);
            indcurvepeak(pt,pk) = pkind1(indcurve1);
            ampcurvepeak(pt,pk) = maxcurve1*sidemult;
        end
    end
end

%now calculate kinematic parameters
tpeak = NaN(size(indpeak));
good = isfinite(indpeak);
tpeak(good) = t(indpeak(good));

tcurvepeak = NaN(size(indcurvepeak));
good = isfinite(indcurvepeak);
tcurvepeak(good) = t(indcurvepeak(good));

smm = nanmean(smm,2);

wavespeed = NaN(1,size(indpeak,2));
for pk = 1:size(indpeak,2)
    p = polyfit(tcurvepeak(2:end,pk),smm(2:end),1);
    wavespeed(pk) = p(1);
end

per = 2*diff(tpeak,[],2);
per(:,size(indpeak,2)) = NaN;
per = nanmean(per);

wavelen = NaN(size(indcurvepeak));
for pk = 1:size(indcurvepeak,2)-1
    for pt = size(indcurvepeak,1)-1:-1:1
        if isfinite(indcurvepeak(pt,pk))
            if (indcurvepeak(pt,pk) >= min(indcurvepeak(:,pk+1))) && ...
                    (indcurvepeak(pt,pk) <= max(indcurvepeak(:,pk+1)))
                %find where on the body the next peak is when this peak is at the tail
                good = isfinite(indcurvepeak(:,pk+1));
                sprevpeak = interp1(indcurvepeak(good,pk+1),smm(good), indcurvepeak(pt,pk));
                wavelen(pt,pk) = 2*(smm(end) - sprevpeak);
            else
                break;
            end
        end
    end
end

amp = NaN(size(indpeak));
amp(:,2:end-1) = (2*amppeak(:,2:end-1) - amppeak(:,1:end-2) - amppeak(:,3:end))/2;

ampcurve = NaN(size(indpeak));
ampcurve(:,2:end-1) = (2*ampcurvepeak(:,2:end-1) - ampcurvepeak(:,1:end-2) - ampcurvepeak(:,3:end))/2;

K.smm = smm;
K.exc = exc;
K.indpeak = indpeak;
K.tpeak = tpeak;
K.sidepeak = sidepeak;
K.indcurvepeak = indcurvepeak;
K.tcurvepeak = tcurvepeak;
K.ampcurvepeak = ampcurvepeak;
K.amp = amp;
K.ampcurve = ampcurve;
K.per = per;
K.wavespeed = wavespeed;
K.wavelen = wavelen;

    
        
    
    
