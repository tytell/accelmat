function acc = process_accel_accel(imu, kin, varargin)

opt.accerr = [0.008 0.008 0.008];

opt = parsevarargin(opt,varargin, 3, 'typecheck',false);

ind = kin.indpeak(end,:);
side = kin.sidepeak(end,:);

good = isfinite(ind);
ind = ind(good);
side = side(good);

tailphase0 = 0.5 * ones(size(ind));
if (side(1) == 'L')
    tailphase0(1) = 0;
end
tailphase0 = cumsum(tailphase0) + 1;
    
tailphase1 = NaN(size(kin.t));

span1 = ind(1):ind(end);
tailphase1(span1) = interp1(kin.t(ind),tailphase0, kin.t(span1),'spline');

span2 = (imu.t >= kin.t(ind(1))) & (imu.t <= kin.t(ind(end)));
tailphase2 = NaN(size(imu.t));
tailphase2(span2) = interp1(kin.t(span1),tailphase1(span1), imu.t(span2));

tpeak = kin.tpeak(end,:);
imudt = imu.t(2) - imu.t(1);
ind2 = round((tpeak - imu.t(1)) / imudt) + 1;

accs = NaN(size(imu.accdyn));
for i = 1:3
    sp = spaps(imu.t(span2), imu.accdyn(span2,i), opt.accerr(i).^2);
    accs(span2,i) = fnval(sp, imu.t(span2));
end

accfwdpk = NaN(1,length(ind2));
accmean = NaN(3,length(ind2));
acciqr = NaN(3,length(ind2));
orientmean = NaN(3,length(ind2));
orientstd = NaN(3,length(ind2));

for i = 2:length(ind2)-1
    if any(~isfinite(ind2(i-1:i+1)))
        continue;
    end
    k = ind2(i-1):ind2(i+1);
    k2 = ind2(i):ind2(i+1);
    
    for j = 1:3
        if any(isfinite(tailphase2(k))) && any(isfinite(accs(k,j)))
            % peak is the 99 perctile forward acceleration for each half
            % tail beat
            if (j == 1)
                pk = prctile(accs(k2,1),99);
                accfwdpk(1,i) = pk;
            end
            % iqr and mean are taken over the whole tail beat
            acciqr(j,i) = iqr(accs(k));            
            accmean(j,i) = nanmean(accs(k,j));
            
            [orientmean(j,i),~,orientstd(j,i)] = angmean(imu.orient(k,j));
        end
    end
end

acc.fwdpk = accfwdpk;
acc.mean = accmean;
acc.iqr = acciqr;
acc.orient = orientmean;
acc.orientstd = orientstd;



