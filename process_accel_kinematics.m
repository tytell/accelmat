function K = process_accel_kinematics(filename, varargin)

opt.view = 'ventral';
opt.mirror = false;     % filmed through a mirror at 45deg (flips y axis)
opt.massperlen = ones(10,1);
opt.smoothdur = 0.5;
opt.fishlen = [];
opt.sampfreq = [];

opt = parsevarargin(opt,varargin, 2);

W = warning;
warning('off','MATLAB:load:variableNotFound');
load(filename,'exs','eys','scale','fishlenmm','haxmmss','haymmss','humms','hvmms',...
    'hxmm','hymm','hxs','hys','t','txmm','tymm','txs','tys','fps');
warning(W);

if exist('exs','var') && ~isempty(exs) && any(~isnan(exs))
    ismidline = true;
else
    ismidline = false;
end
if ~exist('fishlenmm','var')
    if ~isempty(opt.fishlen)
        fishlenmm = opt.fishlen;
    else
        error('No fish length!');
    end
end

good = isfinite(txmm) & isfinite(tymm);

t0 = t(end);
t = t(good) - t0;

if ismidline
    exs = exs(:,good);
    eys = eys(:,good);
end
haxmmss = haxmmss(good);
haymmss = haymmss(good);
humms = humms(good);
hvmms = hvmms(good);
hxmm = hxmm(good);
hymm = hymm(good);
hxs = hxs(good);
hys = hys(good);
txmm = txmm(good);
tymm = tymm(good);
txs = txs(good);
tys = tys(good);

%reorder the points so that they go from head to tail
%first get the head to tail axis
htx = txs - hxs;
hty = tys - hys;
htmag = sqrt(htx.^2 + hty.^2);
htx = htx./htmag;
hty = hty./htmag;

if ismidline
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
else
    mx = [hxs; txs];
    my = [hys; tys];
    npts = 2;
end
nfr = size(mx,2);

if opt.mirror
    my = -my;
end

mxmm = mx*scale;
mymm = my*scale;
dsmm = sqrt(diff(mxmm).^2 + diff(mymm).^2);
smm = [zeros(1,size(mx,2)); cumsum(dsmm)];

if ismidline
    %estimate center of mass position
    s1 = nanmean(smm,2);
    s1 = s1 / s1(end);
    len1 = linspace(0,1,length(opt.massperlen));
    massperlen = interp1(len1,opt.massperlen, s1);

    %integrate to get total mass
    mass = trapz(s1,massperlen);

    %integrate(x(s) * massperlen(s) ds)/area
    % to get COM position
    comx = zeros(1,nfr);
    comy = zeros(1,nfr);
    for i = 1:nfr,
        comx(i) = trapz(s1, mxmm(:,i).*massperlen ./ mass);
        comy(i) = trapz(s1, mymm(:,i).*massperlen ./ mass);
    end;

    %now look for the main body axis
    swimvecx0 = NaN(1,size(mxmm,2));
    swimvecy0 = NaN(1,size(mxmm,2));
    midx = NaN(size(mxmm));
    midy = NaN(size(mymm));

    ind = find(all(isfinite(mxmm)));
    for i = ind
        %straight line interpolant, minimizing the perpendicular distance
        %to the points (doesn't prioritize x or y, and works if we're
        %horizontal or vertical)
        sp1 = spap2(1, 2, smm(:,i)', [mxmm(:,i)-comx(i),mymm(:,i)-comy(i)]');
        xy1 = fnval(sp1,smm([1 end],i)')';

        swimvecx1 = xy1(1,1) - xy1(end,1);
        swimvecy1 = xy1(1,2) - xy1(end,2);
        mag = sqrt(swimvecx1.^2 + swimvecy1.^2);

        swimvecx0(i) = swimvecx1/mag;
        swimvecy0(i) = swimvecy1/mag;
    end

    %smooth over a period longer than a tailbeat
    good = isfinite(swimvecx0);
    swimvecx = NaN(size(swimvecx0));
    swimvecy = NaN(size(swimvecx0));
    swimvecx(good) = get_low_baseline(t(good)',swimvecx0(good)', 1/opt.smoothdur)';
    swimvecy(good) = get_low_baseline(t(good)',swimvecy0(good)', 1/opt.smoothdur)';
    mag = sqrt(swimvecx.^2 + swimvecy.^2);
    swimvecx = swimvecx ./ mag;
    swimvecy = swimvecy ./ mag;

    excx = mxmm - repmat(comx,[npts 1]);
    excy = mymm - repmat(comy,[npts 1]);
else
    comx = hxmm;
    comy = hymm;
    excx = mxmm - repmat(hxmm,[npts 1]);
    excy = mymm - repmat(hymm,[npts 1]);
    
    if nanmedian(abs(txmm - hxmm) ./ abs(tymm - hymm)) > 1
        swimvecx = sign(nanmedian(hxmm - txmm)) * ones(1,nfr);
        swimvecy = zeros(1,nfr);
        ishoriz = true;
    else
        swimvecx = zeros(1,nfr);
        swimvecy = sign(nanmedian(hymm - tymm)) * ones(1,nfr);
        ishoriz = false;
    end
end

switch opt.view
    case 'ventral'
        leftsign = -1;
    case 'dorsal'
        leftsign = 1;
end

exc = -excx.*repmat(swimvecy,[size(excx,1) 1]) + excy.*repmat(swimvecx,[size(excx,1) 1]);

if ismidline 
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

    bad = indpeak == 0;
    indpeak(bad) = NaN;
    amppeak(bad) = NaN;

    good = any(isfinite(indpeak));
    indpeak = indpeak(:,good);
    amppeak = amppeak(:,good);
    sidepeak = sidepeak(:,good);
    npk = size(indpeak,2);

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
else
    [pkamp1,pkind1,pksgn1] = findpeaks2(exc(end,:),'minmax', ...
        'numneighbors',2);

    indpeak = NaN(2,length(pkind1));
    indpeak(end,:) = pkind1;
    amppeak = NaN(size(indpeak));
    amppeak(end,:) = pkamp1;
    npk = length(pkind1);
    sidepeak = spaces(npts,npk);

    sidepeak(end,pksgn1 == 1) = 'L';
    sidepeak(end,pksgn1 == -1) = 'R';
    
    indcurvepeak = NaN(size(indpeak));
    ampcurvepeak = NaN(size(indpeak));
end

%now calculate kinematic parameters
tpeak = NaN(size(indpeak));
good = isfinite(indpeak);
tpeak(good) = t(indpeak(good));

tcurvepeak = NaN(size(indcurvepeak));
good = isfinite(indcurvepeak);
tcurvepeak(good) = t(indcurvepeak(good));

smm = nanmean(smm,2);

%then derivative to get velocity
comvelx = deriv(t,comx);
comvely = deriv(t,comy);

comspeedfwd = comvelx.*swimvecx + comvely.*swimvecy;
comspeedlat = -comvelx.*swimvecy + comvely.*swimvecx;

headdispfwd = [(diff(hxmm).*swimvecx(1:end-1) + diff(hymm).*swimvecy(1:end-1)) NaN];

%now take means for each tail beat
comspeedfwdmn = NaN(1,npk);
comspeedfwdrms = NaN(1,npk);
comspeedlatrms = NaN(1,npk);
headdispfwdmn = NaN(1,npk);
prevheadx = NaN;
prevheady = NaN;
prevswimvecx = NaN;
prevswimvecy = NaN;
for pk = 1:size(indpeak,2)-1
    if isfinite(indpeak(end,pk)) && isfinite(indpeak(end,pk+1))
        btind = indpeak(end,pk)+1:indpeak(end,pk+1);
        
        comspeedfwdmn(pk) = nanmean(comspeedfwd(btind));
        comspeedfwdrms(pk) = rms(comspeedfwd(btind));
        comspeedlatrms(pk) = rms(comspeedlat(btind));
        
        headxmn = nanmean(hxmm(btind));
        headymn = nanmean(hymm(btind));
        
        %displacement of the head from last period to this one, along the
        %last swimming direction
        headdispfwdmn(pk) = (headxmn-prevheadx).*prevswimvecx + ...
            (headymn-prevheady).*prevswimvecy;

        prevswimvecx = nanmean(swimvecx(btind));
        prevswimvecy = nanmean(swimvecy(btind));
        prevheadx = headxmn;
        prevheady = headymn;
    end
end

per = 2*diff(tpeak,[],2);
per(:,size(indpeak,2)) = NaN;
per = nanmean(per);

if ismidline
    wavespeed = NaN(1,size(indpeak,2));
    for pk = 1:size(indpeak,2)
        tc1 = tcurvepeak(2:end,pk);
        s1 = smm(2:end);
        good = isfinite(tc1);
        if sum(good) > 3
            p = polyfit(tc1(good),s1(good),1);
            wavespeed(pk) = p(1);
        end
    end


    wavelen = NaN(size(indcurvepeak));
    for pk = 1:size(indcurvepeak,2)-1
        for pt = size(indcurvepeak,1)-1:-1:1
            if isfinite(indcurvepeak(pt,pk))
                if (indcurvepeak(pt,pk) >= min(indcurvepeak(:,pk+1))) && ...
                        (indcurvepeak(pt,pk) <= max(indcurvepeak(:,pk+1)))
                    %find where on the body the next peak is when this peak is at the tail
                    good = isfinite(indcurvepeak(:,pk+1));
                    if sum(good) >= 2
                        sprevpeak = interp1(indcurvepeak(good,pk+1),smm(good), indcurvepeak(pt,pk));
                        wavelen(pt,pk) = 2*(smm(end) - sprevpeak);
                    end
                else
                    break;
                end
            end
        end
    end
end

amp = NaN(size(indpeak));
amp(:,2:end-1) = (2*amppeak(:,2:end-1) - amppeak(:,1:end-2) - amppeak(:,3:end))/2;

ampcurve = NaN(size(indpeak));
ampcurve(:,2:end-1) = (2*ampcurvepeak(:,2:end-1) - ampcurvepeak(:,1:end-2) - ampcurvepeak(:,3:end))/2;

    
K.t = t;
K.s = smm / fishlenmm;
K.exc = exc / fishlenmm;
K.mx = mxmm / fishlenmm;
K.my = mymm / fishlenmm;
K.exc = exc / fishlenmm;
K.headdispfwd = headdispfwd / fishlenmm;
K.comspeedfwd = comspeedfwd / fishlenmm;
K.comspeedlat = comspeedlat / fishlenmm;
K.indpeak = indpeak;
K.tpeak = tpeak;
K.sidepeak = sidepeak;
K.comspeedfwdmn = comspeedfwdmn / fishlenmm;
K.comspeedfwdrms = comspeedfwdrms / fishlenmm;
K.headdispfwdmn = headdispfwdmn / fishlenmm;
K.amp = amp / fishlenmm;
K.per = per;
if ismidline
    K.curve = curve * fishlenmm;
    K.indcurvepeak = indcurvepeak;
    K.tcurvepeak = tcurvepeak;
    K.ampcurvepeak = ampcurvepeak;
    K.ampcurve = ampcurve * fishlenmm;

    K.wavespeed = wavespeed / fishlenmm;
    K.wavelen = wavelen / fishlenmm;
end

if ~isempty(opt.sampfreq) && (opt.sampfreq ~= fps)
    fprintf('Resampling kinematic data at %f Hz\n', opt.sampfreq);
    
    t0 = t;
    t = t(1):1/opt.sampfreq:t(end);
    
    nold = length(t0);
    nnew = length(t);
    
    Knew.t = t;
    
    fn = fieldnames(K);
    for i = 1:length(fn)
        if ismember(fn{i},{'t'})
            continue;
        end
        if size(K.(fn{i}),2) == nold
            nd = ndims(K.(fn{i}));
            pmt = [2 1 3:nd];
            
            val0 = permute(K.(fn{i}),pmt);
            sz = size(val0);
            n1 = prod(sz(2:end));
            
            val = NaN([nnew sz(2:end)]);
            for j = 1:n1
                good = isfinite(val0(:,j));
                span = (t >= first(t0,good)) & (t <= last(t0,good));
                
                val(span,j) = interp1(t0(good),val0(good,j), t(span));
            end
            
            Knew.(fn{i}) = val;
        else
            Knew.(fn{i}) = K.(fn{i});
        end
    end
    K = Knew;
end
            
    
        
    
    
