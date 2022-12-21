function [C,VP,data] = triangulate_cylinders_wor_opt(PP,ll,Rsol)

nc = size(PP,3);

if nargin <3
%% estimate vp
boffa = zeros(2*nc,4);
for iii = 1:nc
    boffa(2*iii-1,:) = ll(:,1,iii)'*PP(:,:,iii);
    boffa(2*iii,:) = ll(:,2,iii)'*PP(:,:,iii);
end
[~,~,V] = svd(boffa(:,1:3));

Rsol = getRfromV(V);
%r2sol = V(:,end);
%rr = null(r2sol');
%Rsol = [rr(:,1)';r2sol';rr(:,2)'];
%if det(Rsol)<0
   % Rsol = diag([1 1 -1])*Rsol;
%    Rsol = diag([1 -1 1])*Rsol;
%end
end

Tcam = [Rsol' zeros(3,1);0 0 0 1];
Tp = [Rsol zeros(3,1);0 0 0 1];

data.Tcam = Tcam;
data.Tp = Tp;


%% transform to canonical cylinder
PPt = zeros(3,4,nc);

for iii = 1:nc
    PPt(:,:,iii) = PP(:,:,iii)*Tcam;
end

data.PPt = PPt;

%% get 2D-lines in xz-plane
ll2 = zeros(3,2*nc);
for iii = 1:nc
    ltmp = PPt(:,:,iii)'*ll(:,1,iii);
    ll2(:,2*iii-1)=ltmp([1 3 4]);
    ltmp = PPt(:,:,iii)'*ll(:,2,iii);
    ll2(:,2*iii)=ltmp([1 3 4]);    
end


data.ll2 = ll2;



%% find cylinder coeffs

bestdd = solver_triangtransa_opt(ll2);

%disp(bestins)
c0 = [bestdd(1:3) [bestdd(2); bestdd(4:5)] [bestdd(3); bestdd(5); bestdd(6)]];

data.c0 = c0;
if 0
    c0 = -c0/c0(end);
    x0 = c0(1,3);
    r2 = c0(1)+x0^2;
    disp(sqrt(r2));
end

    
C0 = zeros(4);  
C0([1 3 4],[1 3 4]) = c0;

C = Tp'*C0*Tp;
VP = Tcam*[0 1 0 0]';

data.C = C;
data.VP = VP;





