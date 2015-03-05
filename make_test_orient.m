function [imu,corr] = make_test_orient

dist = 0.04;
f = 1;
amp = 20 * pi/180;
angvel0 = 60 * pi/180;
gyrodriftamp = 0.01;
gyronoiseamp = 0.1;
accnoiseamp = 0.1;

t = (0:0.01:4/f)';
n = length(t);

ang = amp*cos(2*pi*f*t);
angvel = [zeros(n,1) -2*pi*f*sin(2*pi*f*t) zeros(n,1)];
angacc = [zeros(n,1) -4*pi^2*f^2*cos(2*pi*f*t) zeros(n,1)];

%ang = t * angvel0;
%angvel = zeros(n,3);
%angvel(:,2) = angvel0 * ones(n,1);
%angacc = zeros(n,3);

r = dist*[cos(ang) zeros(size(ang)) sin(ang)];

velcross = crossProductMatrix(angvel);
acccross = crossProductMatrix(angacc);

acclinw = zeros(n,3);
for i = 1:n
    acclinw(i,:) = r(i,:) * acccross(:,:,i) + (r(i,:) * velcross(:,:,i)) * velcross(:,:,i);
end

rot = zeros(3,3,n);

rot(1,1,:) = cos(ang);
rot(1,3,:) = -sin(ang);
rot(2,2,:) = ones(size(ang));
rot(3,1,:) = sin(ang);
rot(3,3,:) = cos(ang);

g0 = [0 0 -9.8];
gs = zeros(n,3);
acclins = zeros(n,3);
for i = 1:length(t)
    gs(i,:) = g0*rot(:,:,i);
    acclins(i,:) = acclinw(i,:) * rot(:,:,i);
end

accs = acclins + gs;

gyrodrift = cumsum(randn(n,3)*gyrodriftamp);
gyronoise = randn(n,3)*gyronoiseamp;
accnoise = randn(n,3)*accnoiseamp;

imu.t = t;
imu.acc = (accs + accnoise) / 9.8;
imu.gyro = angvel + gyrodrift + gyronoise;

corr.t = t;
corr.acclins = acclins;
corr.angvel = angvel;
corr.gs = gs;
corr.gyrodrift = gyrodrift;
corr.gyronoise = gyronoise;
corr.accnoise = accnoise;

% Outputs the skew symmetric cross product matrix associated with cross
% product operation for a x b => [ax] b where [ax] is a skew symmetric
% matrix. Refer - http://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
% 
% [theMatrix] = crossProductMatrix(theVector)
% Input = theVector - N X 3
% Output = theMatrix - 3 X 3 X N
function C = crossProductMatrix(v)
N = size(v,1);
C = zeros(3, 3, N);
for ii = 1:N
    C(:,:,ii) = [0, -v(ii,3), v(ii,2);
        v(ii,3), 0, -v(ii,1);
        -v(ii,2), v(ii,1), 0];
end
