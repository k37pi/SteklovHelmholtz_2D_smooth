%clear; %close all
set(0,'defaultTextInterpreter','latex'); 
set(0,'defaultAxesFontSize',20)
set(0, 'DefaultLineLineWidth', 2);

a = 1; b = 1; 
l1 = (5); l2 = 2*l1;
l1 = a*l1; l2 = a*l2;
dom_area = l1*l2; perimeter = 2*(l1+l2); 
% dom_area = 1; l1 = 1/3; l2 = 1/l1; 
% perimeter = 2*(l1+l2); 

mu = 2; mu = mu*b/perimeter; 
sc = (ceil(log(mu^2))); 
if sc < 1
    es = 10^sc;  
else
    es = 1e-2;
end    
%% 0 < k = l^2 < mu -----------------
k = 0:es:mu^2; l = sqrt(k);
s = sqrt(mu^2-k);
Y0 = -l.*tan(l*l1/2); Y1 = -s.*tan(s*l2/2);
Y2 = l.*cot(l*l1/2); Y3 = -s.*tan(s*l2/2);
Y4 = -l.*tan(l*l1/2); Y5 = s.*cot(s*l2/2);
Y6 = l.*cot(l*l1/2); Y7 = s.*cot(s*l2/2);

figure()
% subplot(2,2,1)
subplot(3,4,1)
plot(k,Y0,'-'); hold on
plot(k,Y1,'-')
ylim([-1 0])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_1(x,y)$")

% subplot(2,2,2)
subplot(3,4,2)
plot(k,Y2,'-'); hold on
plot(k,Y3,'-')
ylim([-10 10])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_2(x,y)$")

% subplot(2,2,3) 
subplot(3,4,3)
plot(k,Y4,'-'); hold on
plot(k,Y5,'-')
ylim([-10 10])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_3(x,y)$")

%subplot(2,2,4)
subplot(3,4,4)
plot(k,Y6,'-'); hold on
plot(k,Y7,'-')
xlabel("$k$"); ylabel("$\sigma$")
title("$u_4(x,y)$")
legend(["$_x\sigma$","$_y\sigma$"],'Interpreter','latex')

if length(k) > 1
[x1,y1,iout1,jout1] = intersections(k,Y0,k,Y1);
[x2,y2,iout2,jout2] = intersections(k,Y2,k,Y3);
[x3,y3,iout3,jout3] = intersections(k,Y4,k,Y5);
[x4,y4,iout4,jout4] = intersections(k,Y6,k,Y7);
end
%% k = l^2 > mu^2 ----------
if mu == 0 
k = mu^2:es:10; l = sqrt(k);
else
k = mu^2+es:es:10; l = sqrt(k);
end    
s2 = sqrt(k-mu^2);
Y8 = -l.*tan(l*l1/2); Y9 = s2.*tanh(s2*l2/2);
Y10 = l.*cot(l*l1/2); Y11 = s2.*tanh(s2*l2/2);
Y12 = -l.*tan(l*l1/2); Y13 = s2.*coth(s2*l2/2); 
Y14 = l.*cot(l*l1/2); Y15 = s2.*coth(s2*l2/2);

% figure(3)
% subplot(2,2,1)
hpf = subplot(3,4,5);
plot(k,Y8,'-'); hold on
plot(k,Y9,'-')
ylim([min(s2.*tanh(s2*l2/2)) max(s2.*tanh(s2*l2/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_1(x,y)$")
[x5,y5,iout5,jout5] = intersections(k,Y8,k,Y9);

% subplot(2,2,2)
subplot(3,4,6)
plot(k,Y10,'-'); hold on
plot(k,Y11,'-')
ylim([min(s2.*tanh(s2*l2/2)) max(s2.*tanh(s2*l2/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_2(x,y)$")
[x6,y6,iout6,jout6] = intersections(k,Y10,k,Y11);

% subplot(2,2,3)
subplot(3,4,7)
plot(k,Y12,'-'); hold on
plot(k,Y13,'-')
ylim([min(s2.*coth(s2*l2/2)) max(s2.*coth(s2*l2/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_3(x,y)$")
[x7,y7,iout7,jout7] = intersections(k,Y12,k,Y13);

% subplot(2,2,4)
subplot(3,4,8)
plot(k,Y14,'-'); hold on
plot(k,Y15,'-')
ylim([min(s2.*coth(s2*l2/2)) max(s2.*coth(s2*l2/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_4(x,y)$")
legend(["$_x\sigma$","$_y\sigma$"],'Interpreter','latex')
[x8,y8,iout8,jout8] = intersections(k,Y14,k,Y15);

%% 0 > k = -l^2 -----------------
k = -10:es:-es; l = sqrt(-k);
s1 = sqrt(mu^2-k);
Y16 = l.*tanh(l*l1/2); Y17 = -s1.*tan(s1*l2/2);
Y18 = l.*coth(l*l1/2); Y19 = -s1.*tan(s1*l2/2);
Y20 = l.*tanh(l*l1/2); Y21 = s1.*cot(s1*l2/2);
Y22 = l.*coth(l*l1/2); Y23 = s1.*cot(s1*l2/2);

% figure(2)
% subplot(2,2,1)
%d1 = s1*l2/2; 
subplot(3,4,9)
plot(k,Y16,'-'); hold on
plot(k,Y17,'-')
ylim([min(l.*tanh(l*l1/2)) max(l.*tanh(l*l1/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_1(x,y)$")
[x9,y9,iout9,jout9] = intersections(k,Y16,k,Y17);

% subplot(2,2,2)
subplot(3,4,10)
plot(k,Y18,'-'); hold on
plot(k,Y19,'-')
ylim([min(l.*coth(l*l1/2)) max(l.*coth(l*l1/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_2(x,y)$")
[x10,y10,iout10,jout10] = intersections(k,Y18,k,Y19);

%subplot(2,2,3)
subplot(3,4,11)
plot(k,Y20,'-'); hold on
plot(k,Y21,'-')
ylim([min(l.*tanh(l*l1/2)) max(l.*tanh(l*l1/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_3(x,y)$")
[x11,y11,iout11,jout11] = intersections(k,Y20,k,Y21);

% subplot(2,2,4)
subplot(3,4,12)
plot(k,Y22,'-'); hold on
plot(k,Y23,'-')
ylim([min(l.*coth(l*l1/2)) max(l.*coth(l*l1/2))])
xlabel("$k$"); ylabel("$\sigma$")
title("$u_4(x,y)$")
legend(["$_x\sigma$","$_y\sigma$"],'Interpreter','latex')
[x12,y12,iout12,jout12] = intersections(k,Y22,k,Y23);

sgtitle(join(["$\mu = $",mu,", $L_1$=",l1,", $L_2$=",l2,", row 1: $0<k<\mu^2$, row 2: $k>\mu^2$, row 3: $k<0$. Eigen values are at the intersections." ]),fontsize=20)

filename = ['figs/','rectangle_mu',num2str(mu),'_l1',num2str(l1),'_l2',num2str(l2),'.png'];
filename = join(filename,'');
set(findobj(gcf,'type','axes'),'FontSize',20)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
% exportgraphics(gcf,filename,'Resolution',300) 
saveas(gcf,filename);
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% hgexport(ax, strrep(filename,".png",".eps"),...
%     hgexport('factorystyle'), 'Format', 'eps');

if mu ~=0
    ys = [y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12];
    figure()
    plot(ys,'--*'); hold on
    ys = sort(ys);
    plot(ys,'--*')
else
    ys = sort([y5;y6;y7;y8]);
end