function kalman1(t,gyro,acc)

dt = t(2)-t(1);
q1 = 0.1;
q2 = 0.1;
q3 = 0.01;
r1 = 0.1;
r2 = 0.1;

x = zeros(3,length(t));
x(1,1) = atan2(acc(3,1),acc(1,1));

F = [1 dt -dt;
	 0  1   0;
	 0  0   1];
 
Q = [q1  0  0;
	  0 q2  0;
	  0  0 q3];
 
R = [r1  0;
	  0 r2];

H = [1 0 0];

P = 100*eye(3);

for i = 1:length(t)
	x(:,i+1) = F*x(:,i);
 
	P = F*P*F' + Q;
 
	z(1) = atan2(acc(3,i),acc(1,i));
    z(2) = gyro(2,i);
 
	y = z - H*x;
 
	S = H*P*H' + R;
 
	K = P*H*inv(S);
 
	x = x + K*y;
 
	P = (I-K*H)*P;
end 
