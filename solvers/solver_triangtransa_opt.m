function [D,bestmin,DD,bestxyr] = solver_triangtransa_opt(ll)

lille = 1e-8;

ll = ll./sqrt(sum(ll(1:2,:).^2));
la = ll(1,:)';
lb = ll(2,:)';
lc = ll(3,:)';
dd = sum([la.^4, 4*la.^3.*lb, 4*la.^3.*lc, 2*la.^2.*lb.^2, 4*la.^2.*lb.*lc, 2*la.^2.*lc.^2, 4*la.^2.*lb.^2, 8*la.^2.*lb.*lc, ...
        4*la.*lb.^3, 8*la.*lb.^2.*lc, 4*la.*lb.*lc.^2, 4*la.^2.*lc.^2, 4*la.*lb.^2.*lc, 8*la.*lb.*lc.^2, 4*la.*lc.^3, lb.^4, ...
        4*lb.^3.*lc, 2*lb.^2.*lc.^2, 4*lb.^2.*lc.^2, 4*lb.*lc.^3, lc.^4]);
    
sols = solver_triang_opt(dd);
sols(:,sum(abs(imag(sols)))>lille)=[];
a6 = -1;

bestmin = inf;    
nsol = size(sols,2);

DD = [];

for iii = 1:nsol
    a4 = sols(4,iii);
    a5 = sols(5,iii);
    lam1 = sols(6,iii);
    lam2 = sols(7,iii);
    boffa = [2*dd(1) dd(2) dd(3) lam2 - dd(6) + a4*dd(4) + a5*dd(5);...
            dd(2) 2*dd(7) dd(8) lam1 - dd(11) + a4*dd(9) + a5*dd(10);...
            dd(3) dd(8) 2*dd(12) + 2*lam2 a4*dd(13) - dd(15) + a5*dd(14) + a5*lam1;...
            dd(4) dd(9) dd(13) 2*a4*dd(16) - lam2 - dd(18) + a5*dd(17);...
            dd(5) dd(10) dd(14) + lam1 a4*dd(17) - dd(20) + 2*a5*dd(19) - 2*a5*lam2];
     [~,~,V] = svd(boffa);
     solv = V(:,end);
     solv = solv/solv(end);
     a1 = solv(1);
     a2 = solv(2);
     a3 = solv(3);
    
    linvv = [a1^2, a1*a2, a1*a3, a1*a4, a1*a5, a1*a6, a2^2, a2*a3, a2*a4, a2*a5, a2*a6, a3^2, a3*a4, a3*a5, a3*a6, a4^2, a4*a5, a4*a6, a5^2, a5*a6, a6^2];
    
    linvv = linvv/(a1^2);
    
    mini = dd*linvv';
    
    DD = [DD [a1 a2 a3 a4 a5 a6]'];
    
    if mini<bestmin
        bestmin = mini;
        D = [a1 a2 a3 a4 a5 a6]';
    end
end

c0 = [D(1:3) [D(2); D(4:5)] [D(3); D(5); D(6)]];

c0 = -c0/c0(end);
x0 = c0(1,3);
y0 = c0(2,3);
r2 = c0(1)+x0^2;
bestxyr = [x0;y0;sqrt(r2)];

