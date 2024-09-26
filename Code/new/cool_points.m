clear; close all

imname = "flo5";
img = imread(join(["figs/",imname,".jpg"],""));  % Replace with your image file
ims = imshow(img);

RI = imref2d(size(img));
mx = max(RI.XWorldLimits); my = max(RI.YWorldLimits);
RI.XWorldLimits = [0 5];
RI.YWorldLimits = [0 my/mx*5];
hRef = imshow(img,RI);
set(gca,'YDir','reverse');

% figure
% plot(0,0)
% xlim([0 5])
% ylim([0 5])

% place even # of points
[x,y] = getpts;%input("whats the point yo");
x = [x';y']; 
theta = pi;
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
x = rot*x;
a = 0; b = 2*pi;
dt = (b-a)/(length(x)); %cutoff = 0.1; 
t = a:dt:b; % 2N+1 points

x = [x x(:,1)];
x = fliplr(x);

plot(x(1,:),x(2,:),'*'); axis equal

ftx = (x(1,:)); fty = (x(2,:));

fN = 15; 
asx = zeros(1,fN+1); asy = zeros(1,fN+1);
bsx = zeros(1,fN); bsy = zeros(1,fN);
asx(1) = sum(ftx(1:end-1))*dt/(2*pi);
asy(1) = sum(fty(1:end-1))*dt/(2*pi);

for i = 1:fN
    asx(i+1) = sum(ftx(1:end-1).*cos(i*t(1:end-1)))*dt/pi;
    asy(i+1) = sum(fty(1:end-1).*cos(i*t(1:end-1)))*dt/pi;
    bsx(i) = sum(ftx(1:end-1).*sin(i*t(1:end-1)))*dt/pi;
    bsy(i) = sum(fty(1:end-1).*sin(i*t(1:end-1)))*dt/pi;
end    

fsx = @(t) asx(1)+sum(asx(2:end).*cos((1:fN)*t)+bsx(1:end).*sin((1:fN)*t));
fsy = @(t) asy(1)+sum(asy(2:end).*cos((1:fN)*t)+bsy(1:end).*sin((1:fN)*t));
fxr = arrayfun(@(x) fsx(x),t);
fyr = arrayfun(@(x) fsy(x),t);
clf
plot(x(1,:),x(2,:)); hold on;
plot(fxr,fyr); axis equal

writematrix(asx,join(["files/",imname,"_ax.csv"],""))
writematrix(bsx,join(["files/",imname,"_bx.csv"],""))
writematrix(asy,join(["files/",imname,"_ay.csv"],""))
writematrix(bsy,join(["files/",imname,"_by.csv"],""))