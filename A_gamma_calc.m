function A_out = A_gamma_calc(A,t)

[U,S_SVD,V] = svd(A);

Ry = U*V';
Ryz = eye(3);
Ryz(1:2,1:2) = Ry;
q0 = [1,0,0,0];
q1 = Rot2Quat(Ryz);
S = V*S_SVD*V';
alpha = acos(dot(q0,q1))/2;

quat_interp = (q0*sin((1-t)*alpha) + q1*sin(t*alpha))/sin(alpha);
R_interp = Quat2Rot(quat_interp);
R_interp = R_interp(1:2,1:2);
A_out = R_interp*((1-t)*eye(2) + t*S);

end

