function [imu] = load_imu2(filename, varargin)

opt.calib = false;
opt.clickzero = false;
opt.removeweirdblock = true;
opt.resample = true;

% opt = parsevarargin(opt, varargin, 2);

acc = h5read(filename,'/Data/Accel');
gyro = h5read(filename,'/Data/Gyro');
t = h5read(filename,'/Data/t');

istrigger = h5read(filename,'/Data/Trigger');
istrigger = istrigger > 0;
iszero = h5read(filename,'/Data/Zero');
iszero = iszero > 0;

acc_units = h5readatt(filename, '/Data/Accel','Units');
acc_units = char(acc_units');
gyro_units = h5readatt(filename, '/Data/Gyro','Units');
gyro_units = char(gyro_units');
acc_range = h5readatt(filename, '/Data/Accel','full_range');
gyro_range = h5readatt(filename, '/Data/Gyro','full_range');

t_units = h5readatt(filename, '/Data/t','Units');
t_units = char(t_units');

acc_units = regexprep(acc_units, '\\00', '');
gyro_units = regexprep(gyro_units, '\\00', '');
t_units = regexprep(t_units, '\\00', '');

switch t_units
    case {'nanosec','nsec'}
        tscale = 1e9;
    case {'millisec','msec'}
        tscale = 1e3;
    otherwise
        error('Unrecognized time unit: %s\n', t_units);
end

t = t/tscale;

good = t ~= 0;
if (opt.removeweirdblock)
    t1 = t(good);
    skipind = find(diff(t1) > 1);
    if ~isempty(skipind)
        t0 = t1(skipind(end)+1);
        good = good & (t >= t0);
    end
end

% % % keyboard
% % % 
% % % t = t(good);
% % % acc = acc(good,:);
% % % gyro = gyro(good,:);
% % % istrigger = istrigger(good,:);
% % % iszero = iszero(good,:);
% % % 
% % % if any(istrigger)
% % %     t0 = first(t,istrigger);
% % % else
% % %     t0 = first(t,isfinite(t));
% % % end
% % % t = t-t0;
% % % 
% % % if opt.calib
% % %     if opt.clickzero
% % %         fig = figure;
% % %         mb = msgbox('Click three zero times','Click','help','non-modal');
% % %         plot(t,acc);
% % %         vertplot(t(iszero),'k--');
% % %         
% % %         [tzero,~] = ginputb(3);
% % %         iszero = false(size(t));
% % %         for i = 1:3
% % %             [~,ind] = min(abs(t-tzero(i)));
% % %             iszero(ind) = true;
% % %         end
% % %         close(mb);
% % %         close(fig);
% % %     end
% % %     
% % %     axord = inputdlg({'First zero:', 'Second zero:', 'Third zero:'}, ...
% % %         'Order of axes',1,{'Y','-Z','X'});
% % %     axord = regexp(axord, '([+-]?)([XYZxyz])', 'tokens', 'once');
% % %     
% % %     axsign = ones(3,1);
% % %     axname = char(zeros(3,1));
% % %     for i = 1:3
% % %         if ~isempty(axord{i}{1}) && (axord{i}{1} == '-')
% % %             axsign(i) = -1;
% % %         end
% % %         axname(i) = upper(axord{i}{2});
% % %     end
% % %     
% % %     [~,ord] = sort(axname);
% % %     axsign = axsign(ord);
% % %     
% % %     gax = acc(iszero,:)';
% % %     gax = gax(:,ord);
% % %     gax = gax .* repmat(axsign', [3 1]);
% % %     
% % %     basis = gramschmidt(gax);
% % %     %check for right handed-ness
% % %     %the Z axis should be equal to the cross product of the X and Y axes.
% % %     %because of small numerical issues, it's sometimes not exactly equal,
% % %     %so we check that they're in the same direction
% % %     assert(sum(cross(basis(:,1),basis(:,2)) .* basis(:,3)) > 0.9);
% % %     
% % %     calib.chip2world = basis;
% % %     calib.world2chip = basis';
% % % end
% % % 
% % % if opt.resample
% % %     newsamp = diground(1./nanmean(diff(t)), 10);
% % %     fprintf('Resampling at %d Hz\n', newsamp);
% % %     t0 = t;
% % %     t = (t0(1):1/newsamp:t0(end))';
% % %     acc = interp1(t0,acc, t);
% % %     gyro = interp1(t0,gyro, t);
% % % end

imu.t = t;
imu.acc = acc;% * calib.chip2world;
imu.acc_units = acc_units;
imu.acc_range = acc_range;
imu.gyro = gyro;% * calib.chip2world;
imu.gyro_units = gyro_units;
imu.gyro_range = gyro_range;


