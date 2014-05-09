function imu = load_imu(filename, varargin)

opt.removeweirdblock = true;
opt.extrapoints = 20;
opt = parsevarargin(opt, varargin, 2);

acc = h5read(filename, '/Data/Accel');
gyro = h5read(filename, '/Data/Gyro');
t = h5read(filename, '/Data/t');
istrigger = h5read(filename, '/Data/Trigger');
istrigger = istrigger > 0;
iszero = h5read(filename, '/Data/Zero');
iszero = iszero > 0;

acc_range = h5readatt(filename, '/Data/Accel','full_range');
acc_units = h5readatt(filename, '/Data/Accel','Units');
gyro_range = h5readatt(filename, '/Data/Gyro','full_range');
gyro_units = h5readatt(filename, '/Data/Gyro','Units');
t_units = h5readatt(filename, '/Data/t','Units');

acc_units = char(acc_units');
acc_units = regexprep(acc_units,'\\00','');
gyro_units = char(gyro_units');
gyro_units = regexprep(gyro_units,'\\00','');
t_units = char(t_units');
t_units = regexprep(t_units,'\\00','');

switch t_units
    case {'nsec','nanosec'}
        tscale = 10^9;
    case {'microsec','usec'}
        tscale = 10^6;
    case {'millisec','msec'}
        tscale = 10^3;
    otherwise
        error('Unknown units for time');
end

t = t/tscale;
good = t ~= 0;

if (opt.removeweirdblock)
    t1 = t(good);
    isskip = diff(t1) > 1;
    
    if (any(isskip))
        tweird = t1(last(isskip)+1);
        
        good = good & (t >= tweird);
    end
end

ind = find(iszero);
expts = opt.extrapoints + 1;
acc0 = zeros(expts,3,length(ind));
gyro0 = zeros(expts,3,length(ind));
off = (1:expts)-floor(expts/2);

for i = 1:length(ind)
    acc0(:,:,i) = acc(ind(i)+off,:);
    gyro0(:,:,i) = acc(ind(i)+off,:);
end

ttrig = first(t,istrigger);

t = t(good);
acc = acc(good,:);
gyro = gyro(good,:);

if (any(ttrig))
    t = t - ttrig;
else
    t = t - t(1);
end

imu = struct('t',t,'acc',acc,'acc_units',acc_units,'acc_range',acc_range, ...
    'gyro',gyro, 'gyro_units',gyro_units, 'gyro_range',gyro_range, ...
    'ttrig',ttrig, ...
    'acc0', acc0, 'gyro0',gyro0, 't0',t(ind));


