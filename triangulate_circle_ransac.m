function [D,bestxyr,inliers,bestins] = triangulate_circle_ransac(ll,bnd,niters)
% function [D,bestxyr,inliers,bestins] = triangulate_circle_ransac(ll,bnd,niters)
%
% input
% ll: 3xn 2d-lines 
% bnd: inliers bnd
% niters: nr of iterations
% output 
% D: 6x1 vector of best dual conic
% bestxyr: 3x1 vector containing best position (x,y) and radius 
% inliers: 1xn boolean vector of inliers
% bestins: nr of found inliers

if nargin<2
    bnd = 1e-5;
end
if nargin<3
    niters = 1000;
end

n = size(ll,2);

alldata = [ll(1,:).^2; 2*ll(1,:).*ll(2,:); 2*ll(1,:).*ll(3,:);...
           ll(2,:).^2; 2*ll(2,:).*ll(3,:); ll(3,:).^2];
alldata = alldata./sqrt(sum(alldata.^2));
alldata = alldata';

alldata1 = alldata(1:2:end,:);
alldata2 = alldata(2:2:end,:);


bestins = 0;
% Ransac

for iii = 1:niters
    ids = randperm(n,3);
    %dd = solver_triangtransa2(ll(:,ids(1)),ll(:,ids(2)),ll(:,ids(3)));
    dd = triangulate_circle_minimal(ll(:,ids(1)),ll(:,ids(2)),ll(:,ids(3)));
    
    %insi = sum(abs(alldata*dd)<bnd);
    insi = sum((abs(alldata1*dd)+abs(alldata2*dd))<bnd);
    
    [maxinsi,idsi] = max(insi);
    if maxinsi>bestins
        bestins = maxinsi;
        bestdd = dd(:,idsi);
    end
end

inliers = abs(alldata*bestdd)<bnd;
   
D = bestdd;
c0 = [bestdd(1:3) [bestdd(2); bestdd(4:5)] [bestdd(3); bestdd(5); bestdd(6)]];
c0 = -c0/c0(end);
x0 = -c0(1,3);
y0 = -c0(2,3);
r = sqrt(c0(1)+x0^2);
bestxyr = [x0;y0;r];

