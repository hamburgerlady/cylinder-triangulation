function D = triangulate_circle_minimal(l1,l2,l3)
% function D = triangulate_circle_minimal(l1,l2,l3)
%
% input
% l1,l2,l3 3x1  vectors of 3 2d-lines
% output
% D: 6xnsol dual conic solutions


lille = 1e-8;

% constr x3*x5-x2*x6, x3^2-x5^2-x1*x6+x4*x6
% constr r=1 x5^2-x4x6-x6^2 x3x5-x2x6 x3x4-x2x5+x3x6 x3^2-x1x6-x6^2 x2x3-x1x5-x5x6 x2^2-x1x4-x1x6-x4x6-x6^2


Asol = null([l1(1)^2 2*l1(1)*l1(2) 2*l1(1)*l1(3), l1(2)^2, 2*l1(2)*l1(3), l1(3)^2;...
    l2(1)^2 2*l2(1)*l2(2) 2*l2(1)*l2(3), l2(2)^2, 2*l2(2)*l2(3), l2(3)^2;...
    l3(1)^2 2*l3(1)*l3(2) 2*l3(1)*l3(3), l3(2)^2, 2*l3(2)*l3(3), l3(3)^2]);
    


sol = solver_triang_null(Asol);
sol(:,sum(abs(imag(sol)))>lille)=[];

D = Asol*[sol;ones(1,size(sol,2))];
