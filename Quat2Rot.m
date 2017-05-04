function R = Quat2Rot(q)
R = zeros(3);
sq = q.*q;
invs = 1.0/sum(sq);
R(1,1) = ( sq(2) - sq(3) - sq(4) + sq(1)) *invs;
R(2,2) = (-sq(2) + sq(3) - sq(4) + sq(1)) *invs;
R(3,3) = (-sq(2) - sq(3) + sq(4) + sq(1)) *invs;

tmp1 = q(2)*q(3);
tmp2 = q(4)*q(1);
R(2,1) = 2.0 * (tmp1 + tmp2)*invs ;
R(1,2) = 2.0 * (tmp1 - tmp2)*invs ;

tmp1 = q(2)*q(4);
tmp2 = q(3)*q(1);
R(3,1) = 2.0 * (tmp1 - tmp2)*invs ;
R(1,3) = 2.0 * (tmp1 + tmp2)*invs ;
tmp1 = q(3)*q(4);
tmp2 = q(2)*q(1);
R(3,2) = 2.0 * (tmp1 + tmp2)*invs ;
R(2,3) = 2.0 * (tmp1 - tmp2)*invs ;    
end

