function out = trerp(p,q, t)
%rigid triangle interpolation
q_interp = q;

M = zeros(6,6);
for i=1:3
    M(2*i-1,:) = [p(i,1),p(i,2),0,0,1,0];
    M(2*i,:) = [0,0,p(i,1),p(i,2),0,1];
end

q_prime = q';
q_prime = q_prime(:);

Al = ((M'*M)^-1)*M'*q_prime;

A = [Al(1),Al(2);Al(3),Al(4)];

[U,S_SVD,V] = svd(A);

Ry = U*V;
Ryz = eye(3);
Ryz(1:2,1:2) = Ry;
q0 = Rot2Quat(eye(3));
q1 = Rot2Quat(Ryz);
S = V'*S_SVD*V;

alpha = acos(Ry(1))/2;
quat_interp = (q0*sin((1-t)*alpha) + q1*sin(t*alpha))/sin(alpha);
R_interp = Quat2Rot(quat_interp);
R_interp = R_interp(1:2,1:2);
newA = R_interp*((1-t)*eye(2) + t*S);

for j = 1:3
    q_interp(j,:) = newA*p(j,:)' + t*Al(5:6);
end
out = q_interp;
end

