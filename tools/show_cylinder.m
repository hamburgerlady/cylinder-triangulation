function ptc = myshow_cylinder3(D,VP,PP,imseq,params)



if nargin<5
    params.hmax = 20;
    params.hmin = 0;
    params.nnh = 15;
    params.nnth = 15;
    params.color = 'g';
    params.linewidth = 1;
else
    if ~isfield(params,'hmax')
        params.hmax = 70;
    end
    if ~isfield(params,'hmin')
        params.hmin = 0;
    end
    if ~isfield(params,'nnh')
        params.nnh = 15;
    end
    if ~isfield(params,'nnth')
        params.nnth = 15;
    end
    if ~isfield(params,'color')
        params.color = 'g';
    end
    if ~isfield(params,'linewidth')
        params.linewidth = 1;
    end
end


    
nnh = params.nnh;
hmax = params.hmax;
hmin = params.hmin;
r2 = VP(1:3);
R = [r2 null(r2')];
R = R(:,[2 1 3]);
if det(R)<0
    R(:,3) = -R(:,3);
end

T = [R zeros(3,1);0 0 0 1];
D2 = T'*D*T;
d2 = D2([1 3 4],[1 3 4]);
c2 = inv(d2);

x0 = pflat(null(c2(1:2,:)));
T2 = [eye(2) x0(1:2);0 0 1];
c2t = T2'*c2*T2;
c2t = -c2t/c2t(end);
r = (sqrt(1/c2t(1)));

th = linspace(0,2*pi,500);
hh = linspace(hmin,hmax,nnh);

XX = [];
XXd = [];
hold on
for iii = 1:nnh
    X = R*([r*cos(th);0*th+hh(iii);r*sin(th)]+[x0(1);0;x0(2)]);
    Xd = R*([-cos(th);0*th;-sin(th)]);
    XX = [XX X];
    XXd = [XXd Xd];
end
nnp = 100000;
idx = randperm(size(XX,2),nnp);
XX = XX(:,idx);

XXe = [XX;ones(1,nnp)];

XXd = XXd(:,idx);

[M,N,~] = size(imseq(:,:,:,1));


P1 = squeeze(PP(1,:,:));
P2 = squeeze(PP(2,:,:));
P3 = squeeze(PP(3,:,:));

p1 = P1'*XXe;
p2 = P2'*XXe;
p3 = P3'*XXe;

dsc = (p1./p3-N/2).^2+(p2./p3-M/2).^2;


nim = size(PP,3);

C0 = zeros(3,nim);
d0 = zeros(3,nim);
for iii = 1:nim
    Pi = PP(:,:,iii);
    [~,Pin] = rq(Pi);
    d0(:,iii) = Pin(3,1:3)';
    tmp = pflat(null(Pi));
    C0(:,iii) = tmp(1:3);
end
    


dsc0 = d0'*XXd;
okids = dsc0>0.5;

dsc(~okids)=inf;


[MM,idi] = min(dsc);

idi(MM>((M/2)^2+(N/2)^2))=0;

col = uint8(zeros(size(XX)));
for iii = 1:nim
   idii = find(idi==iii);
   im = imseq(:,:,:,iii);
   P = PP(:,:,iii);
   Xi = XX(:,idii);
   u = P*[Xi;ones(1,size(Xi,2))];
    u = round(pflat(u));
    imr = im(:,:,1);
    img = im(:,:,2);
    imb = im(:,:,3);
    [M,N] = size(imr);
    u(2,u(2,:)>M)=nan;
    u(2,u(2,:)<1)=nan;
    u(1,u(1,:)>N)=nan;
    u(1,u(1,:)<1)=nan;

    ind = sub2ind([M N],u(2,:),u(1,:));
    ind(isnan(ind))=1;
    col(:,idii) = [imr(ind);img(ind);imb(ind)];
end

ptc = pointCloud(XX','Color',col');

