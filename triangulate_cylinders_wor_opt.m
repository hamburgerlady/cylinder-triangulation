function [C,VP,data] = triangulate_cylinders_wor_opt(PP,ll)


% rectify
[ll2,data] = rectify_cylinders(PP,ll);

% find cylinder coeffs

[Dopt,bestmin,DD,bestxyr_opt] = triangulate_circle_opt(ll2);

c0 = [Dopt(1:3) [Dopt(2); Dopt(4:5)] [Dopt(3); Dopt(5); Dopt(6)]];
data.c0 = c0;
data.x0 = bestxyr_opt(1);
data.y0 = bestxyr_opt(2);
data.r = bestxyr_opt(3);
data.bestmin = bestmin;
data.allsols = DD;



C0 = zeros(4);  
C0([1 3 4],[1 3 4]) = c0;

% transform to original coordinate system
Tp = data.Tp;
Tcam = data.Tcam;

C = Tp'*C0*Tp;
VP = Tcam*[0 1 0 0]';

data.C = C;
data.VP = VP;





