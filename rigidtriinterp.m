p = [1,2;0.5,4;2,4];
q = [2,1;2,6;4,5];
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

figure(1);
clf;
%axis([0, 6, 0, 6]);
t = 1/(24);
alpha = acos(Ry(1))/2;

for i = 0:t:1
    quat_interp = (q0*sin((1-i)*alpha) + q1*sin(i*alpha))/sin(alpha);
    R_interp = Quat2Rot(quat_interp);
    R_interp = R_interp(1:2,1:2);
    newA = R_interp*((1-i)*eye(2) + i*S);
    for j = 1:3
        q_interp(j,:) = newA*p(j,:)' + i*Al(5:6);
    end
    axis([0, 6, 0, 6]);
    hold on;
    DrawBetweenTwoPoints(q_interp(1,:),q_interp(2,:));
    DrawBetweenTwoPoints(q_interp(2,:),q_interp(3,:));
    DrawBetweenTwoPoints(q_interp(3,:),q_interp(1,:));
    drawnow;
end