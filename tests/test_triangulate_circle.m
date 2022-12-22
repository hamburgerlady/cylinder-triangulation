%% generate data
addpath ../matlab/tools/

% cylinder
x0 = rand(2,1)*10;
r = 1+rand;

% tangent lines
nlines = 20;
th = rand(1,nlines)*2*pi;
N0 = [cos(th);sin(th)];
P0 = x0 + r*N0;
L0 = [N0;-P0(1,:).*N0(1,:)-P0(2,:).*N0(2,:)];



bruce = 0.005;
outs = 5;
ll = L0+randn(size(L0))*bruce;
ll(:,randperm(nlines,outs))=randn(3,outs);




%% test least squares
[Dopt,bestmin,DD,bestxyr_opt] = triangulate_circle_opt(ll);
disp([bestxyr_opt [x0;r]])

%% test minimal
ll = ll./sqrt(sum(ll(1:2,:).^2));
[Dmin,bestxyr_min,inliers,bestins] = triangulate_circle_ransac(ll,1e-3,1000);
disp(bestins)

%% test least squares on inlier set
[Dopt2,bestmin2,DD2,bestxyr_opt2] = triangulate_circle_opt(ll(:,inliers));
disp([bestxyr_opt2 [x0;r]])


%% plots
figure(1); % least squares
clf
hold on
rital(L0);

th0 = linspace(0,2*pi,100);
lh = plot(x0(1)+r*cos(th0),x0(2)+r*sin(th0));
set(lh,'LineWidth',2);

lh = plot(bestxyr_opt(1)+bestxyr_opt(3)*cos(th0),bestxyr_opt(2)+bestxyr_opt(3)*sin(th0),'--');
set(lh,'LineWidth',2);

%plot(P0(1,:),P0(2,:),'*');

axis equal
axis([-5 15 -5 15])


figure(2); % least squares
clf
hold on
rital(L0);

th0 = linspace(0,2*pi,100);
lh = plot(x0(1)+r*cos(th0),x0(2)+r*sin(th0));
set(lh,'LineWidth',2);

lh = plot(bestxyr_min(1)+bestxyr_min(3)*cos(th0),bestxyr_min(2)+bestxyr_min(3)*sin(th0),'--');
set(lh,'LineWidth',2);

%plot(P0(1,:),P0(2,:),'*');

axis equal
axis([-5 15 -5 15])

figure(3); % least squares on inliers
clf
hold on
rital(L0);

th0 = linspace(0,2*pi,100);
lh = plot(x0(1)+r*cos(th0),x0(2)+r*sin(th0));
set(lh,'LineWidth',2);

lh = plot(bestxyr_opt2(1)+bestxyr_opt2(3)*cos(th0),bestxyr_opt2(2)+bestxyr_opt2(3)*sin(th0),'--');
set(lh,'LineWidth',2);

%plot(P0(1,:),P0(2,:),'*');

axis equal
axis([-5 15 -5 15])

