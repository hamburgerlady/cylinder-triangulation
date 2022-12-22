function [ll2,data] = rectify_cylinders(PP,ll)

nc = size(PP,3);

boffa = zeros(2*nc,4);
for iii = 1:nc
    boffa(2*iii-1,:) = ll(:,1,iii)'*PP(:,:,iii);
    boffa(2*iii,:) = ll(:,2,iii)'*PP(:,:,iii);
end
[~,~,V] = svd(boffa(:,1:3));

Rsol = getRfromV(V);


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

