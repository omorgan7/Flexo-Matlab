function [P_inv,A] = trerp(p,q)
%rigid triangle interpolation
q_interp = q;

M = zeros(6,6);
for i=1:3
    M(2*i-1,:) = [p(i,1),p(i,2),1,0,0,0];
    M(2*i,:) = [0,0,0,p(i,1),p(i,2),1];
end

q_prime = q';
q_prime = q_prime(:);
P_inv = (M'*M)^-1*M';
Al = P_inv*q_prime;
P_inv = [P_inv(1:2,:);P_inv(4:5,:)];
A = [Al(1),Al(2);Al(4),Al(5)];
end

