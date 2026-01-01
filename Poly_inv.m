function [p_inv1, out1, p_inv2, out2, out2_check] = Poly_inv(p,in)

p_inv1=0;
out1 = 0;
Ce = p;
% p1 + p2*X + p3*X^2 + p4*X^3

% Convert the polynomial to a companion matrix
%p_inv1 = [1 Ce(1) Ce(1)^2 Ce(1)^3;0 Ce(2) 2*Ce(1)*Ce(2) 3*Ce(1)*Ce(2);0 Ce(3) 2*Ce(1)*Ce(3)+Ce(2)^2 3*Ce(1)*(Ce(1)*Ce(3)+Ce(2)^2);0 Ce(4) 2*Ce(1)*Ce(4)+2*Ce(2)*Ce(3) 3*Ce(1)^2*Ce(4)+6*Ce(1)*Ce(2)*Ce(3)+Ce(2)^3]\[0;1;0;0];

p_inv2(1) = -p(1)/p(2) - p(3)*p(1)^2/p(2)^3 + p(4)*p(1)^3/p(2)^4;
p_inv2(2) = 1/p(2) + 2*p(1)*p(3)/p(2)^3 - 3*p(1)^2*p(4)/p(2)^4;
p_inv2(3) = -p(3)/p(2)^3 + 3*p(1)*p(4)/p(2)^4;
p_inv2(4) = -p(4)/p(2)^4;
% test
% T = [-100 -10 -1 0 1 10 100];
% 
% Out1 = p(1) + p(2)*T + p(3)*T.^2 + p(4)*T.^3
% 
% out1 = p_inv1(1) + p_inv1(2)*in + p_inv1(3)*in.^2 + p_inv1(4)*in.^3;

 out2 = p_inv2(1) + p_inv2(2)*in + p_inv2(3)*in.^2 + p_inv2(4)*in.^3;

 out2_check = ((in-p(1))/p(2)) - p(3)/p(2)*((in-p(1))/p(2)).^2 - p(4)/p(2)*((in-p(1))/p(2)).^3;



end