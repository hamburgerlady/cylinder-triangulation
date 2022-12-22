%% get lines
bp = 'tests/';
dd = dir([bp '*jpg']);
nim = length(dd);
load([bp 'kungshuset_data.mat']);



%% plot extracted lines 


for iii = 1:nim
    im = imread([bp dd(iii).name]);
    if iii == 1
        [M,N,NC] = size(im);
        imseq = uint8(zeros(M,N,NC,nim));
    end
    imseq(:,:,:,iii) = im;
    hold off
    imshow(im)
    hold on
    rital(ll(:,:,iii));
    title(num2str(iii));
    drawnow
    pause
end


%% triangulate

[C,VP,data] = triangulate_cylinders_wor_opt(PPn,lln);


%% show results 
figure(1);
clf
hold on
plotconic(m2v(data.c0),'b');
for iii = 1:nim
    tmp = pflat(null(data.PPt(:,:,iii)));
    plot(tmp(1),tmp(3),'*');
end

axis equal
axis([-10 10 -5 10])


figure(2);
clf
pt0 = pointCloud(XX');
pcshow(pt0);
hold on
params.hmin = -7;
params.hmax = 0.5;
params.nnh = 1000;
params.linewidth = 2;
ptc = show_cylinder(C,VP,PP,imseq,params);
pcshow(ptc)



    
    


