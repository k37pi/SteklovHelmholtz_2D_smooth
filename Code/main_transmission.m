clear; close all; 
%% inputs -----------------------------------------------
a = 0; b = 2*pi; 
Ns = [240];%,420]; % N = half the number of intervals, set more than mu^2
curves_list = ["" + ...
    "Circle"; "Ellipse"; "kite";"coolcurve";"disk_def";
    "acorn";"star";"rectangle";"egg";"opt_disk";"butterfly";"cavity";"bunimovich"];%;"peanut"];
% cool curve = (a1 cos (b1 t) + a2 sin (b2 t))ang2a[cost, sint] 

% circle = 1, ellipse = 2,  kite = 3, coolcurve = 4, 
% disk_def = 5, acorn = 6 , star = 7, square = 8, egg = 9
% opt_disk = 10, butterfly = 11, cavity = 12, bunimovich = 13
curve_name = curves_list(9); curve_name = lower(curve_name);

% curve params
disk_r = 7; 
ellipse_a = 1; ellipse_b = 1;
cool_ab = [2,2,3];
square_a = 1; square_b = 2; epps_egg = 0.2;
%cool_ab = randi([1,70],2,1); 
%cool_ab = [cool_ab;randi([1,100],1,1)]; 
% cool_ab = [cool_ab;randi([cool_ab(1),100],1,1)]; 
% c1 = (cool_ab(3)^2-cool_ab(1)^2)/(4*pi^2*cool_ab(1)^2);
% b1 = sqrt((1+sqrt(1+4*c1))/2);
% cool_ab(2) = randi([floor(0.5*(cool_ab(1)+cool_ab(3))),100],1,1);

%cool_ab = [33;3;68];

diskdef_k = 3.7188;%-5.30538947734613e-06; 
len1 = 0;

% inside grid
m = 150;

RI = 2; % refractive index for transmission
%mu_min = 9.45; mu_max = 9.45;
MU = 3.5175;
mu_min = 10; mu_max = 10.5;
mus = mu_min:(mu_max-mu_min)/299:mu_max;
step = diff(mus); step = step(1);
%mus = 3.38239;
%mus = 5.9997;4.9026;

%min_mu = 0; counter = 1; % # digits after decimal
count = 0; total_count = length(mus)*length(Ns); % progress

%while(numel(num2str(min_mu)) < 10)

s_mins = zeros(length(mus),1); min_inds = zeros(length(mus),1);

for ns = 1:length(Ns)
N = Ns(ns); all_S = zeros(2*N,length(mus));
%--------------------------------------------------------
%% curve points and parametrization ---------------------
dt = (b-a)/(2*N); %cutoff = 0.1; 
t = a:dt:b; % 2N+1 points

switch curve_name
    case "circle"
        [x,dx,d2x,~,len] = curve(curve_name,t,len1,disk_r);
        curve_param = join([curve_name,disk_r],"");
    case "ellipse"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,ellipse_a,ellipse_b);
        curve_param = join([curve_name,num2str(ellipse_a),'\_',num2str(ellipse_b)],"");
    case "coolcurve"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,cool_ab);
        %curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3)),'\_',num2str(cool_ab(4)),'\_'],"");
        curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3))],"");
    case "disk_def"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,diskdef_k);  
        curve_param = join([curve_name,diskdef_k],"");
    case "rectangle"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,square_a,square_b);  
        curve_param = join([curve_name,num2str(square_a),'\_',num2str(square_b)],"");
    case "egg"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,epps_egg);  
        curve_param = join([curve_name,num2str(epps_egg)],"");
    otherwise
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1);
        curve_param = curve_name;
end

fig_curvename = strrep(curve_param,".","dot");

xi = linspace(min(x(1,:)),max(x(1,:)),m) ;
yi = linspace(min(x(2,:)),max(x(2,:)),m) ;
[X,Y] = meshgrid(xi,yi) ;
idx = inpolygon(X(:),Y(:),x(1,:),x(2,:)) ;
points_inside = transpose([X(idx),Y(idx)]);

%points_inside_close = 0.95*x(:,1:2:end-1);

points_inside_close = x(:,1:2:end-1) - 0.05*nx(:,1:2:end-1);

points_inside = [points_inside,points_inside_close];

[Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
[Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

%figure(); 
plot(x(1,:),x(2,:),'LineWidth',1.4);   hold on; 
quiver(x(1,37),x(2,37),nx(1,37),nx(2,37),'LineWidth',1.4); hold on;
plot(points_inside_close(1,:),points_inside_close(2,:),'o');
% text(0,0,'\Omega')
% legend(["$M = \partial \Omega$","$\nu = $ outward unit normal"],'Interpreter','latex')
% ylim([-1.5 2])
axis equal

phis = zeros(2*N,length(mus));

for ms = 1:length(mus)

for ri = 1:2
if ri == 1
Mu = mus(ms);%5.07;%212; % hopefully mu^2 is not a D/N ev 
else
Mu = sqrt(RI)*Mu;
end    

if Mu > 0 
    mu = Mu;
else
    mu = -1i*Mu; % or +i Mu ?
end    

[Tau,T] = meshgrid(t);
[Xtau,Xt] = meshgrid(x(1,:));
[Ytau,Yt] = meshgrid(x(2,:));
[dXtau,dXt] = meshgrid(dx(1,:));
[dYtau,dYt] = meshgrid(dx(2,:)); 
[d2Xtau,d2Xt] = meshgrid(d2x(1,:));
[d2Ytau,d2Yt] = meshgrid(d2x(2,:));
%----------------------------------------

%% R, X for all pairs t, tau
tic
X = dYt.*(Xtau-Xt)-dXt.*(Ytau-Yt);
R = sqrt((Xt-Xtau).^2+(Yt-Ytau).^2);
XR_ratio = X./R;
toc
% ----------------------------------------

%% L, M ---------------------------------------
tic
dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);
dpnorm_ratio = dp2norm./dp1norm;
logsin_term = log(4*sin((T-Tau)/2).^2);
cL = (1i*mu/2)*dpnorm_ratio;
cL1 = -mu/(2*pi)*dpnorm_ratio;
L = cL.*XR_ratio.*besselh(1,mu*R); % note that besselh(1,0) blows up so when R = 0 it's useless.

L1 = cL1.*XR_ratio.*besselj(1,mu*R);
% if Mu < 0
%     L1 = 1i*L1;
% end    

L2 = L-L1.*logsin_term;
if curve_name == "rectangle"
    L(isnan(L)) = 0; 
    L2(isnan(L2)) = 0; 
else
    L(isnan(L)) = diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2)); 
    L2(isnan(L2)) = diag(L); 
end
    
L1(isnan(L1)) = 0;


M = (1i/2)*dp2norm.*besselh(0,mu*R);
M1 = -1/(2*pi)*dp2norm.*besselj(0,mu*R);
% if Mu < 0
%     M1 = 1i*M1;
% end    
% if Mu < 0
%     M = (1i/2)*dp2norm.*besselk(0,abs(mu)*R);
%     %M = (1i/2)*dp2norm.*(besseli(0,mu*R)+1i*besselk(0,mu*R));
%     M1 = -1/(2*pi)*dp2norm.*besselj(0,abs(mu)*R);    
% end    
M2 = M-M1.*logsin_term;
mc = 1i/2-double(eulergamma)/pi;

% if Mu > 0
if curve_name == "rectangle"
    M2(isnan(M2)) = 0;
else    
    M2(isnan(M2)) = diag(dp1norm.*(mc-log(mu/2*dp1norm)/pi));
end
    % else
%     M2(isinf(M2)) = diag(dp1norm.*(mc-log(abs(mu)/2*dp1norm)/pi));
% end    
toc
%-------------------------------------------

%% Trig. poly. -----------------------------
% f^tp(t,tau) = 
% (a_0/2,a_1,...,a_{N-1},a_N/2).(1,cos(tau),...,cos([N-1]tau),cos([N]tau))+
% (b_1,...,b_{N-1}).(sin(tau),...,sin([N-1]tau)) =
% a_vec.cos_0N_vec + b_vec.sin_0N_vec
t_scaled = (0:2*N).*T; % the first column is 0*t, second is 1*t, so on.
cos_vecs = cos(t_scaled); % l is fixed in columns, t_j is fixed in rows.
sin_vecs = sin(t_scaled);

Tau_sub_0N_scaled = t_scaled(:,1:N+1); 
cos_0N_vecs = cos_vecs(:,1:N+1);
sin_0N_vecs = sin_vecs(:,2:N);
tic

%L1 tp
L1_tp = ftp(L1,N,cos_vecs,sin_vecs);

% L2 tp
L2_tp = ftp(L2,N,cos_vecs,sin_vecs);

%M1 tp
M1_tp = ftp(M1,N,cos_vecs,sin_vecs);

% M2 tp
M2_tp = ftp(M2,N,cos_vecs,sin_vecs);
toc
% ------------------------------------------

%% RN weights, same for any integrand (L1,M1 for us) ---------------
tic
diff_t_tau = T-Tau;

% if N < 501
%     rep_diff = repmat(diff_t_tau,[1 1 N-1]); % create 3d matrix with copies of diff
%     rep_diff = rep_diff.*reshape(1:N-1,1,1,[]); % multiply each matrix with next integer
%     rep_diff_cos = cos(rep_diff); 
%     rep_diff_cos_scaled = rep_diff_cos./reshape(1:N-1,1,1,[]);
%     RN_mat = -2*pi*sum(rep_diff_cos_scaled,3)/N - pi*cos(N*(diff_t_tau))/N^2;
% else
    % RN_mat = zeros(size(diff_t_tau));
    % ncount = ceil(N/500);   
    % rep_diff1 = repmat(diff_t_tau,[1 1 500]); % create 3d matrix with copies of diff
    % 
    % %tic
    % for ncn = 1:(ncount-1)
    %     inds = (500*(ncn-1)+1):(500*ncn);
    %     %rep_diff = repmat(diff_t_tau,[1 1 300]); % create 3d matrix with copies of diff
    %     % rep_diff = rep_diff1.*reshape(inds,1,1,[]); % multiply each matrix with next integer
    %     % rep_diff_cos = cos(rep_diff); 
    %     % rep_diff_cos_scaled = rep_diff_cos./reshape(inds,1,1,[]);
    %     rep_diff_cos_scaled = cos(rep_diff1.*reshape(inds,1,1,[]))./reshape(inds,1,1,[]);
    %     RN_mat = RN_mat -2*pi*sum(rep_diff_cos_scaled,3)/N - pi*cos(N*(diff_t_tau))/N^2;
    % end 
    % %toc
    % inds = (500*(ncount-1)+1):(N-1);
    % rep_diff1 = rep_diff1(:,:,1:length(inds)); % create 3d matrix with copies of diff
    % % rep_diff = rep_diff.*reshape(inds,1,1,[]); % multiply each matrix with next integer
    % % rep_diff_cos = cos(rep_diff); 
    % % rep_diff_cos_scaled = rep_diff_cos./reshape(inds,1,1,[]);
    % rep_diff_cos_scaled = cos(rep_diff1.*reshape(inds,1,1,[]))./reshape(inds,1,1,[]);
    % RN_mat = RN_mat -2*pi*sum(rep_diff_cos_scaled,3)/N - pi*cos(N*(diff_t_tau))/N^2;


    RN_mat = cos(diff_t_tau);
    for k = 2:(N-1)
        RN_mat = RN_mat+cos(k*diff_t_tau)/k; 
    end    
    RN_mat = -2*pi*RN_mat/N - pi*cos(N*(diff_t_tau))/N^2;
% end
toc

%-------------------------------------------

%% Matrices ------------------------------
A = RN_mat.*L1_tp+pi/N*L2_tp; % (A + I) is the effect of K'
A(:,end) = []; A(end,:) = [];

I = diag(ones(2*N,1));

B = RN_mat.*M1_tp+pi/N*M2_tp; % B is the effect of the S 
B(:,end) = []; B(end,:) = [];


%% Transmission --------------------------------
if ri == 1
    N_mu1 = (I+A)*(B\I);
    B1 = B; A1 = A;
else
    N_mu2 = (I+A)*(B\I);
    N_mu = N_mu1 - N_mu2;
end
end

tic; [U,S,V] = svd(N_mu); toc
S = diag(S); [s_min,min_ind] = min(S); 
all_S(:,ms) = S;
% Use V column corresponding to smallest sv as the "null" vector
% this is the density we will use to compute the solution
fix = min_ind;
phi_min = V(:,fix);
%norm(N_mu*phi_min) % figure(); plot(vecnorm(N_mu*V,2,1))

Xindiff = Xin-Xtauin; Yindiff = Yin-Ytauin;
Rin = sqrt(Xindiff.^2+Yindiff.^2);

% HankRin = besselh(0,mu*Rin)-besselh(0,mu*sqrt(RI)*Rin); % rows correspond to fixed point inside 
% Recon_integrand = HankRin.*transpose(phi_min).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
% recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
% recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
% recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
% recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
% recint_sin = recint_sin1*sin_vecs(:,2:N)';
% recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
% solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
%solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 

HankRin1 = besselh(0,mu*Rin); % rows correspond to fixed point inside 
Recon_integrand = HankRin1.*transpose(phi_min).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
%Recon_integrand = HankRin1.*transpose((B1\I)*phi_min).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
solution_fixed_point1 = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
%solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 

HankRin2 = besselh(0,mu*sqrt(RI)*Rin); % rows correspond to fixed point inside 
Recon_integrand = HankRin2.*transpose(phi_min).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
%Recon_integrand = HankRin2.*transpose((B\I)*phi_min).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
solution_fixed_point2 = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
%solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 

xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));
%[thetalin,rholin] = cart2pol(xlin,ylin);

%[Xin2,Yin2] = meshgrid(xlin,ylin);
%[Rr, theta] = meshgrid(thetalin,rholin);
% Uin = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
% RealUin = real(Uin);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
% %RealUin = RealUin/max(abs(solution_fixed_point));
% ImagUin = imag(Uin);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');
% %ImagUin = ImagUin/max(abs(solution_fixed_point));
% bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
% Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;


[Xin2,Yin2] = meshgrid(xlin,ylin);
%[Rr, theta] = meshgrid(thetalin,rholin);
Uin1 = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point1,Xin2,Yin2,'cubic');
RealUin1 = real(Uin1);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
%RealUin = RealUin/max(abs(solution_fixed_point));
ImagUin1 = imag(Uin1);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');
%ImagUin = ImagUin/max(abs(solution_fixed_point));
bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;
%figure()
ax1=subplot(3,2,[1 3]); surf(Xin2,Yin2,RealUin1)
shading interp; view(0,90); colorbar; colormap(ax1,viridis)
title("$n_{RI}$ = $1$",'Interpreter','latex')


[Xin2,Yin2] = meshgrid(xlin,ylin);
%[Rr, theta] = meshgrid(thetalin,rholin);
Uin2 = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point2,Xin2,Yin2,'cubic');
RealUin2 = real(Uin2);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
%RealUin = RealUin/max(abs(solution_fixed_point));
ImagUin2 = imag(Uin2);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');
%ImagUin = ImagUin/max(abs(solution_fixed_point));
bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;

ax2=subplot(3,2,[2 4]); surf(Xin2,Yin2,RealUin2)
shading interp; view(0,90); colorbar; colormap(ax2,magma)
%axis equal
title(join(["$n_{RI}$ = ",num2str(RI)],""),'Interpreter','latex')

subplot(3,2,5:6);plot(t(1:end-1),phi_min,"LineWidth",2)
xlabel("$t$",'Interpreter','latex')
ylabel("$\phi$",'Interpreter','latex')
xlim([0 2*pi])

sgtitle(join(["Real part of transmission eigenfunctions, $\mu_T$ = ",num2str(mus(ms))],""),'interpreter','latex','FontSize',20)

set(findobj(gcf,'type','axes'),'FontSize',20)
filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,curve_param,"evecs_",mus(ms),"_N",num2str(N)],'');
filename = strrep(filename,'.',"dot");
exportgraphics(gcf,join([filename,".png"],""),'Resolution',300) 
saveas(gcf,filename);

% 
% subplot(1,3,3); surf(Xin2,Yin2,RealUin1-RealUin2)
% shading interp; view(0,90); colorbar
% 
% figure();
% subplot(1,3,1); plot(real(N_mu*phi_min))
% subplot(1,3,2); plot(real(N_mu1*phi_min))
% subplot(1,3,3); plot(real(N_mu2*phi_min))
% 
% figure()
% subplot(1,2,1); plot(t(1:end-1),real(phi_min))
% subplot(1,2,2); surf(Xin2,Yin2,RealUin1-RealUin2)
% shading interp; view(0,90); colorbar
%figure()
subplot(1,2,1)
plot(t(1:end-1),B1*phi_min)
subplot(1,2,2)
plot(t(1:end-1),B*phi_min)

%figure()
%plot(abs(phi_min)); title("$|\phi_{\min}|$",'Interpreter','latex')

% figure(); plot3(x(1,1:end-1),x(2,1:end-1),real(phi_min))

wantinds = size(recint,1)-size(points_inside_close,2)+1;
% figure(); plot3(points_inside(1,wantinds:end),points_inside(2,wantinds:end),real(solution_fixed_point(wantinds:end)),'--o')
% title("Real $\phi_\Gamma$",'Interpreter','latex')

Recon_integrand = (HankRin1.*transpose((B1\I)*phi_min)-HankRin2.*transpose((B\I)*phi_min)).*dp2norm(min_ind,(1:end-1));
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+recint_sin; recint = recint(:,1:end-1);

B_aug = 1i*pi/(4*N)*recint(wantinds:end,:); 
%sum(B_aug*phi_min-solution_fixed_point(wantinds:end) < 1e-10) == N

C_aug = [N_mu;B_aug];
[Q,R] = qr(C_aug); QA = Q(1:size(N_mu,1),:);
QB = Q((size(N_mu,1)+1:end),:);

tic; [U2,S2,V2] = svd(QA,"econ"); toc
S2 = diag(S2); [s_min2,min_ind2] = min(S2);

norm(QA*R*phi_min) % (s_min) should be ~ 0 vector

s_mins(ms) = s_min;
phis(:,ms) = phi_min;
min_inds(ms) = min_ind;

count = count+1;
disp(round(count/total_count*100,2))
%disp(round(ms/length(mus)*100,2))
end

%figure(); 
plot(mus,s_mins,'--o','LineWidth',2)
xlabel('$\mu$','Interpreter','latex')
ylabel('$s_{\min}$','Interpreter','latex')
title(join(["Minimum singular values of $C(\mu)$ for ",curve_param,", $n_{RI}$= ",RI,", N = ",num2str(N)],""),'Interpreter','latex')
filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,fig_curvename,"min_sig_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".png"],"");
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
set(findobj(gcf,'type','axes'),'FontSize',18)
% exportgraphics(gcf,filename,'Resolution',300)  
% savefig(gcf,strrep(filename,".png",".fig"))

filename = 'files/transmission/';%['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,fig_curvename,"_ri_",RI,"min_phi_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".csv"],"");
writematrix(phis,filename)

% fucking save the min sig with mus as csv tf is wrong with you ? so much
filename = 'files/transmission/';%['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,fig_curvename,"_ri_",RI,"minsigmu_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".csv"],"");
writematrix([mus',s_mins],filename)

filename = 'files/transmission/';%['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,fig_curvename,"_ri_",RI,"min_inds_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".csv"],"");
writematrix(min_inds,filename)
end

%min_mu = mus(s_mins == min(s_min));

%end


