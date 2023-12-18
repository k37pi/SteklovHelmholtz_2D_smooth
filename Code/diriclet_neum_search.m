tic
clear; close all; 
%% inputs -----------------------------------------------
a = 0; b = 2*pi; 
%mus = 2.350:0.001:2.404;%[1, 5, 7, 10]; % hopefully mu^2 is not a D/N ev 
%mus = 2.412:-0.0001:2.4049;
%mus = 3.82:0.0001:3.8316;
%mus = 3.831639:0.000001:3.831699;
%mus = 7.015539:0.000001:7.015599;
%mus = 11.6199:-1e-5:11.6197;
%mus = 11.7914:1e-5:11.7916;
%mus = 9.7609:1e-5:9.7611;
%mus = 2.935:1e-4:2.945;
%mus = 1.0194:1e-5:1.0197;
%mus = 4.624:1e-5:4.6249;
%mus = 3.3003:1e-5:3.3005;
mus = 0.1:0.5*1e-2:2;3.4226:1e-5:3.4228;2.4134:1e-5:2.4136;2:1e-1:4;%1.96695:1e-5:1.96705;3.5:1e-2:4;%0.1:1e-2:2;%5.063:1e-4:5.065;%1.68:1e-3:1.7;%2:1e-2:5;%1.2024:1e-5:1.2026;    

%mus = 3.84:-0.0001:3.8318;
N = 200; % N = half the number of intervals, set more than mu^2
curves_list = ["" + ...
    "Circle"; "Ellipse"; "kite";"coolcurve";"disk_def";
    "acorn";"star";"rectangle";"egg";"opt_disk";"butterfly";"cavity";"bunimovich"];%;"peanut"];
% cool curve = (a1 cos (b1 t) + a2 sin (b2 t))ang2a[cost, sint] 

% circle = 1, ellipse = 2,  kite = 3, coolcurve = 4, 
% disk_def = 5, acorn = 6 , star = 7, square = 8, egg = 9
% opt_disk = 10, butterfly = 11, cavity = 12, bunimovich = 13
curve_name = curves_list(1); curve_name = lower(curve_name);
% curve params
ellipse_a = cosh(2); ellipse_b = sinh(2); 
disk_r = 3;5;3/(2*pi);
len1 = 0;
sig = 2;5.5235; %0.886;%-0.5751;%-1;

%Mu = 9.7610;%11.7915;%11.6198;%7.0156;%3.8317;%2.4048;

% inside grid
m = 150;
fix = 1;
% evecs_in = cell(length(mus),1); evecs_M = cell(length(mus),1);
% evs_all = cell(length(mus),1);
% Xin_all = cell(length(mus),1); Yin_all = cell(length(mus),1);
Smins = zeros(length(mus),1);
Smink = zeros(length(mus),1);
Sminr = zeros(length(mus),1);
Sminr1 = zeros(length(mus),1);


%--------------------------------------------------------
for ms = 1:(length(mus))

mu = mus(ms);
%% curve points and parametrization ---------------------
dt = (b-a)/(2*N); 
t = a:dt:b; % 2N+1 points

switch curve_name
    case "circle"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,disk_r);
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

dom_area = polyarea(x(1,:)',x(2,:)');

xi = linspace(min(x(1,:)),max(x(1,:)),m) ;
yi = linspace(min(x(2,:)),max(x(2,:)),m) ;
[X,Y] = meshgrid(xi,yi) ;
idx = inpolygon(X(:),Y(:),x(1,:),x(2,:)) ;
points_inside = transpose([X(idx),Y(idx)]);

[Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
[Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

% figure(); plot(x(1,:),x(2,:)); hold on; plot(points_inside(1,:),points_inside(2,:),'o');
% axis equal

[Tau,T] = meshgrid(t);
[Xtau,Xt] = meshgrid(x(1,:));
[Ytau,Yt] = meshgrid(x(2,:));
[dXtau,dXt] = meshgrid(dx(1,:));
[dYtau,dYt] = meshgrid(dx(2,:)); 
[d2Xtau,d2Xt] = meshgrid(d2x(1,:));
[d2Ytau,d2Yt] = meshgrid(d2x(2,:));
%----------------------------------------

%% R, X for all pairs t, tau
%tic
X = dYt.*(Xtau-Xt)-dXt.*(Ytau-Yt);
R = sqrt((Xt-Xtau).^2+(Yt-Ytau).^2);
XR_ratio = X./R;
%toc
% ----------------------------------------

%% L, M ---------------------------------------
%tic
dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);
dpnorm_ratio = dp2norm./dp1norm;
logsin_term = log(4*sin((T-Tau)/2).^2);

cL = (1i*mu/2)*dpnorm_ratio;
cL1 = -mu/(2*pi)*dpnorm_ratio;
L = cL.*XR_ratio.*besselh(1,mu*R); % note that besselh(1,0) blows up so when R = 0 it's useless.
L1 = cL1.*XR_ratio.*besselj(1,mu*R);
L2 = L-L1.*logsin_term;
L(isnan(L)) = diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2)); 
L1(isnan(L1)) = 0;
L2(isnan(L2)) = diag(L); 

M = (1i/2)*dp2norm.*besselh(0,mu*R);
M1 = -1/(2*pi)*dp2norm.*besselj(0,mu*R);
M2 = M-M1.*logsin_term;
mc = 1i/2-double(eulergamma)/pi;
M2(isnan(M2)) = diag(dp1norm.*(mc-log(mu/2*dp1norm)/pi));
%toc
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
%tic

%L1 tp
L1_tp = ftp(L1,N,cos_vecs,sin_vecs);

% L2 tp
L2_tp = ftp(L2,N,cos_vecs,sin_vecs);

%M1 tp
M1_tp = ftp(M1,N,cos_vecs,sin_vecs);

% M2 tp
M2_tp = ftp(M2,N,cos_vecs,sin_vecs);
%toc
% ------------------------------------------

%% RN weights, same for any integrand (L1,M1 for us) ---------------
%tic
diff_t_tau = T-Tau;

if N < 501
    rep_diff = repmat(diff_t_tau,[1 1 N-1]); % create 3d matrix with copies of diff
    rep_diff = rep_diff.*reshape(1:N-1,1,1,[]); % multiply each matrix with next integer
    rep_diff_cos = cos(rep_diff); 
    rep_diff_cos_scaled = rep_diff_cos./reshape(1:N-1,1,1,[]);
    RN_mat = -2*pi*sum(rep_diff_cos_scaled,3)/N - pi*cos(N*(diff_t_tau))/N^2;
else
    RN_mat = cos(diff_t_tau);
    for k = 2:(N-1)
        RN_mat = RN_mat+cos(k*diff_t_tau)/k; 
    end    
    RN_mat = -2*pi*RN_mat/N - pi*cos(N*(diff_t_tau))/N^2;
end
%toc

%-------------------------------------------

%% Matrices ------------------------------
A = RN_mat.*L1_tp+pi/N*L2_tp; % (A + I) is the effect of K'
A(:,end) = []; A(end,:) = [];

I = diag(ones(2*N,1));

B = RN_mat.*M1_tp+pi/N*M2_tp; % B is the effect of the S 
B(:,end) = []; B(end,:) = [];


%% Dirichlet eigen search --------------------------------
% [U,S1,V] = svd(B); toc
% S = diag(S1); [s_min,min_ind] = min(S); 
% Smins(ms) = s_min; 
% 
% [U,S2,V] = svd(A+I); toc
% S = diag(S2); [s_min,min_ind] = min(S); 
% Smink(ms) = s_min; 

%[U,S,V] = svd(A+I-sig*B); toc
%S = diag(S); [s_min,min_ind] = min(S); 
%Sminr(ms) = s_min; 

ee = eigs(A+I-sig*B,1,"smallestabs"); 
Sminr(ms) = ee; 
%Sminr1(ms) = min([diag(S1);diag(S2)]); 

    %all_S(:,ms) = S;
    % phi_min = V(:,min_ind); evec_dir = B/2*phi_min;
    % HankRin = besselh(0,Mu*Rin); % rows correspond to fixed point inside 
    % Recon_integrand = HankRin.*transpose(phi_min).*dp2norm(min_ind,(1:end-1)); % make this a tp and then sum and scale
    % recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
    % recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2; 
    % recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
    % recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
    % recint_sin = recint_sin1*sin_vecs(:,2:N)';
    % recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
    % solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
    % %solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 
    % xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
    % ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));
    % %[thetalin,rholin] = cart2pol(xlin,ylin);
    % 
    % [Xin2,Yin2] = meshgrid(xlin,ylin);
    % %[Rr, theta] = meshgrid(thetalin,rholin);
    % 
    % Uin = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
    % RealUin = real(Uin);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
    % RealUin = RealUin/max(abs(solution_fixed_point));
    % ImagUin = imag(Uin);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');
    % ImagUin = ImagUin/max(abs(solution_fixed_point));
    % bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
    % Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;
    disp(ms/length(mus)*100)
end    

%% plots --------
figure()
if curve_name == "circle"
Nn = -20:20;
X1=[]; X2=[];X3=[];
for nns=1:length(Nn)
n=Nn(nns);
fun = @(x) x.*(besselj(n-1,disk_r*x)-besselj(n+1,disk_r*x))/2-sig*besselj(n,disk_r*x);
fun1 = @(x) besselj(n,disk_r*x);
fun2 = @(x) (besselj(n-1,disk_r*x)-besselj(n+1,disk_r*x))/2;

%fun = @(x) (besselj(n-1,x)-besselj(n+1,x))/2;

lb = min(mus);             % Set a lower bound for the function.
ub = max(mus);     % Set an upper bound for the function.
nn = 1000;
x = NaN*ones(nn,1);             % Initializes x.
x1 = NaN*ones(nn,1);             % Initializes x.
x2 = NaN*ones(nn,1);             % Initializes x.
starting_points=linspace(lb,ub,nn);
for i=1:nn
        % Look for the zeros in the function's current window.
        x(i)=fzero(fun, starting_points(i));
        x1(i)=fzero(fun1, starting_points(i));
        x2(i)=fzero(fun2, starting_points(i));
end
X1=[X1; uniquetol(x)];
X2=[X2; uniquetol(x1)];
X3=[X3; uniquetol(x2)];
end
X1 = sort(uniquetol(X1)); X2 = sort(uniquetol(X2)); X3 = sort(uniquetol(X3));

% figure()
% plot(lb:1e-3:ub,fun(lb:1e-3:ub))
% yline(0)
% 
end

subplot(3,1,1)
plot(mus,Smins,'--*','LineWidth',1.5)
xlabel("$\mu$",'Interpreter','latex')
ylabel("$s^D_{min}$",'Interpreter','latex')
%title(join(["Possible Laplace Dirichlet eigenvalues of ",curve_param],""),'Interpreter','latex')
title(join(["Dirichlet-Laplace"],""),'Interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',25)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
if curve_name == "circle"
%xline([2.4048	3.8317	5.1356	6.3802	7.5883	8.7715...
%5.5201	7.0156	8.4172	9.7610	11.0647...
%8.6537	10.1735	11.6198 11.7915	9.9361 11.0864])
xline(X2);%/disk_r)
xlim([min(mus) max(mus)])
xticks(X2)
xtickangle(75)
end

subplot(3,1,2)
plot(mus,Smink,'--*r','LineWidth',1.5)
xlabel("$\mu$",'Interpreter','latex')
ylabel("$s^N_{min}$",'Interpreter','latex')
%title(join(["Possible Laplace Neumann eigenvalues of ",curve_param],""),'Interpreter','latex')
title(join(["Neumann-Laplace"],""),'Interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',25)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
if curve_name == "circle"
% xline([3.8317	1.8412	3.0542	4.2012	5.3175	6.4156 ...
%     7.0156 7.5013 5.3314	6.7061	8.0152	9.6474 9.2824	10.5199 ...
%     10.1735	8.5363	9.9695	11.3459 10.7114 11.7060 11.7349 ])
xline(X3);%/disk_r)
xlim([min(mus) max(mus)])
xticks(X3)
xtickangle(75)
end

subplot(3,1,3)
plot(mus,Sminr,'--*','color',"#7E2F8E",'LineWidth',1.5)
xlabel("$\mu$",'Interpreter','latex')
ylabel("$s^R_{min}$",'Interpreter','latex')
%title(join(["Possible Laplace Neumann eigenvalues of ",curve_param],""),'Interpreter','latex')
title(join(["Robin-Laplace, $\sigma$ = ",sig],""),'Interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',25)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
if curve_name == "circle"
xline(X1);%/disk_r)
xlim([min(mus) max(mus)])
xticks(X1)
xtickangle(75)
end
% 4.3887

% xline([3.8317	1.8412	3.0542	4.2012	5.3175	6.4156 ...
%     7.0156 7.5013 5.3314	6.7061	8.0152	9.6474 9.2824	10.5199 ...
%     10.1735	8.5363	9.9695	11.3459 10.7114 11.7060 11.7349 ])

% if curve_name == "circle"
% sgtitle(join([curve_param,", vertical lines correspond to true eigenvalues $\mu^2$"],""),'interpreter','latex','fontsize',20)
% else
% sgtitle(join([curve_param," eigenvalues are $\mu^2$"],""),'interpreter','latex','fontsize',20)
% end    

filename = join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/LapDNRev_minmu_',num2str(min(mus)),'maxmu_',num2str(max(mus)),curve_param,'1.png'],"");
filename = strrep(filename,"\","");
% exportgraphics(gcf,filename,'Resolution',300)
% savefig(strrep(filename,".png",".fig"))

% filename = "/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/lapneumdirev_minmu_1maxmu_5butterfly.fig";
% savefig(filename)
% 
% filename = "/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/lapneumdirev_minmu_1maxmu_5butterfly.png";
% exportgraphics(gcf,filename,'Resolution',300)

% mu_sminr = [Sminr,mus'];
% 
% I = 1:length(mus); % all lam points
% J = I(2:end-1); % interior points
% S = Sminr; 
% J = J( S(J)<S(J-1) & S(J)<S(J+1) ); % local minima
% mus(J)

%% Sign flip, RBF interpolation and zeros
figure()
mu_min = min(mus);
mu_max = max(mus);

% Dirichlet
I = 1:length(mus); % all lam points
J = I(2:end-1); % interior points
S = Smins'; 
J = J( S(J)<S(J-1) & S(J)<S(J+1) ); % local minima
J= J+ (S(J-1)>S(J+1)); % points where sign changes
K = 0*I; K(J) = 1;
S = S.*(-1).^cumsum(K); % introduce sign flips
subplot(3,1,1)
%plot(mus,S,'*',"linewidth",2), hold on % plot signed angle function
%plot([min(mus) max(mus)],[0 0],'-k') % plot lambda axis

rbf_points = mus;
rbf_f = S;

[Xtau_rbf,Xt_rbf] = meshgrid(rbf_points(1,:));
R_rbf = abs(Xt_rbf-Xtau_rbf);

rbf_c = 0.05;
%rbf = @(t) t; drbf = @(t) 1;
%rbf = @(t) t.^3; drbf = @(t) 3*t; % dividing with t for the x/r(t) part. 
%rbf = @(t) t.^2.*log(t+rbf_c); drbf = @(t) 2*log(t+rbf_c)+t./(t+rbf_c);
rbf = @(t) 1./sqrt(t.^2+rbf_c^2); drbf = @(t) -1./(t.^2+rbf_c^2).^1.5;
%rbf = @(t) 1./(t.^2+rbf_c^2); drbf = @(t) -1./(t.^2+rbf_c^2).^2;

% rbf = @(x,y) 1./(x.^2+rbf_c^2).*(1./(y.^2+rbf_c^2)); 
% drbfx = @(x,y) (-2*x)./(x.^2+rbf_c^2).^2*1./(y.^2+rbf_c^2);
% drbfy = @(x,y) (-2*y)./(y.^2+rbf_c^2).^2*1./(x.^2+rbf_c^2);

%rbf = @(t) exp(-rbf_c*t.^2); drbf = @(t) -2*rbf_c*exp(-rbf_c*t.^2);

bf = 0;
if bf == 1
    A_rbf = rbf((Xt_rbf-Xtau_rbf),(Yt_rbf-Ytau_rbf));
else
    A_rbf = rbf(R_rbf);
end
A_rbf = (A_rbf+A_rbf')/2;

rbf_mat = A_rbf;

%rbf_rhs = [rbf_f;zeros(p_basis_size,1)];
rbf_rhs = rbf_f';
%rbf_coeffs = rbf_mat\rbf_lhs;
rbf_coeffs = lsqminnorm(rbf_mat,rbf_rhs);
%poly_coefs = rbf_coeffs((end-p_basis_size+1):end);

points = mu_min:(mu_max-mu_min)/9999:mu_max;

% rbf_fun = zeros(length(points),1);
% for i = 1:length(points)
%     point = points(i); 
%     for j = 1:length(rbf_points)
%         distr = norm(point-rbf_points(j)); 
%         rbf_fun(i) = rbf_fun(i)+rbf_coeffs(j)*rbf(distr);
%     end    
% end

rbf_fun = RBF(points,rbf,rbf_coeffs,rbf_points);
fun = @(x) RBF(x,rbf,rbf_coeffs,rbf_points);

lb = mu_min;             % Set a lower bound for the function.
ub = mu_max;          % Set an upper bound for the function.
nn = 1000;
x = NaN*ones(nn,1);             % Initializes x.
starting_points=linspace(lb,ub,nn);
for i=1:nn
        % Look for the zeros in the function's current window.
        x(i)=fzero(fun, starting_points(i));
end
mu_D = uniquetol(x);
plot(points,rbf_fun,'LineWidth',3)
xline(mu_D,'color',"k","linewidth",2);
yline(0)
xticks(mu_D)
xtickangle(75)
xlim([mu_min mu_max])

xlabel("$\mu$","Interpreter","latex")
ylabel("RBF","Interpreter","latex")
title("Dirichlet-Laplace",'Interpreter','latex')

% Neumann
I = 1:length(mus); % all lam points
J = I(2:end-1); % interior points
S = Smink'; 
J = J( S(J)<S(J-1) & S(J)<S(J+1) ); % local minima
J= J+ (S(J-1)>S(J+1)); % points where sign changes
K = 0*I; K(J) = 1;
S = S.*(-1).^cumsum(K); % introduce sign flips
%hold off, 
subplot(3,1,2)
%plot(mus,S,'*r',"linewidth",2), hold on % plot signed angle function
%plot([min(mus) max(mus)],[0 0],'-k') % plot lambda axis

rbf_points = mus;
rbf_f = S;

[Xtau_rbf,Xt_rbf] = meshgrid(rbf_points(1,:));
R_rbf = abs(Xt_rbf-Xtau_rbf);

rbf_c = 0.05;
%rbf = @(t) t; drbf = @(t) 1;
%rbf = @(t) t.^3; drbf = @(t) 3*t; % dividing with t for the x/r(t) part. 
%rbf = @(t) t.^2.*log(t+rbf_c); drbf = @(t) 2*log(t+rbf_c)+t./(t+rbf_c);
rbf = @(t) 1./sqrt(t.^2+rbf_c^2); drbf = @(t) -1./(t.^2+rbf_c^2).^1.5;
%rbf = @(t) 1./(t.^2+rbf_c^2); drbf = @(t) -1./(t.^2+rbf_c^2).^2;

% rbf = @(x,y) 1./(x.^2+rbf_c^2).*(1./(y.^2+rbf_c^2)); 
% drbfx = @(x,y) (-2*x)./(x.^2+rbf_c^2).^2*1./(y.^2+rbf_c^2);
% drbfy = @(x,y) (-2*y)./(y.^2+rbf_c^2).^2*1./(x.^2+rbf_c^2);

%rbf = @(t) exp(-rbf_c*t.^2); drbf = @(t) -2*rbf_c*exp(-rbf_c*t.^2);

bf = 0;
if bf == 1
    A_rbf = rbf((Xt_rbf-Xtau_rbf),(Yt_rbf-Ytau_rbf));
else
    A_rbf = rbf(R_rbf);
end
A_rbf = (A_rbf+A_rbf')/2;

rbf_mat = A_rbf;

%rbf_rhs = [rbf_f;zeros(p_basis_size,1)];
rbf_rhs = rbf_f';
%rbf_coeffs = rbf_mat\rbf_lhs;
rbf_coeffs = lsqminnorm(rbf_mat,rbf_rhs);
%poly_coefs = rbf_coeffs((end-p_basis_size+1):end);

points = mu_min:(mu_max-mu_min)/9999:mu_max;

% rbf_fun = zeros(length(points),1);
% for i = 1:length(points)
%     point = points(i); 
%     for j = 1:length(rbf_points)
%         distr = norm(point-rbf_points(j)); 
%         rbf_fun(i) = rbf_fun(i)+rbf_coeffs(j)*rbf(distr);
%     end    
% end

rbf_fun = RBF(points,rbf,rbf_coeffs,rbf_points);
fun = @(x) RBF(x,rbf,rbf_coeffs,rbf_points);

lb = mu_min;             % Set a lower bound for the function.
ub = mu_max;          % Set an upper bound for the function.
nn = 1000;
x = NaN*ones(nn,1);             % Initializes x.
starting_points=linspace(lb,ub,nn);
for i=1:nn
        % Look for the zeros in the function's current window.
        x(i)=fzero(fun, starting_points(i));
end
mu_N = uniquetol(x);
plot(points,rbf_fun,'color',"r",'LineWidth',3)
xline(mu_N,'color',"k","linewidth",2);
yline(0)
xticks(mu_N)
xtickangle(75)
xlim([mu_min mu_max])

xlabel("$\mu$","Interpreter","latex")
ylabel("RBF","Interpreter","latex")
title("Neumann-Laplace",'Interpreter','latex')

% Robin
I = 1:length(mus); % all lam points
J = I(2:end-1); % interior points
S = Sminr'; 
J = J( S(J)<S(J-1) & S(J)<S(J+1) ); % local minima
J= J+ (S(J-1)>S(J+1)); % points where sign changes
K = 0*I; K(J) = 1;
S = S.*(-1).^cumsum(K); % introduce sign flips
subplot(3,1,3)
%hold off, plot(mus,S,'*','color',"#7E2F8E","linewidth",2), hold on % plot signed angle function
%plot([min(mus) max(mus)],[0 0],'-k') % plot lambda axis

rbf_points = mus;
rbf_f = S;

[Xtau_rbf,Xt_rbf] = meshgrid(rbf_points(1,:));
R_rbf = abs(Xt_rbf-Xtau_rbf);

rbf_c = 0.05;
%rbf = @(t) t; drbf = @(t) 1;
%rbf = @(t) t.^3; drbf = @(t) 3*t; % dividing with t for the x/r(t) part. 
%rbf = @(t) t.^2.*log(t+rbf_c); drbf = @(t) 2*log(t+rbf_c)+t./(t+rbf_c);
rbf = @(t) 1./sqrt(t.^2+rbf_c^2); drbf = @(t) -1./(t.^2+rbf_c^2).^1.5;
%rbf = @(t) 1./(t.^2+rbf_c^2); drbf = @(t) -1./(t.^2+rbf_c^2).^2;

% rbf = @(x,y) 1./(x.^2+rbf_c^2).*(1./(y.^2+rbf_c^2)); 
% drbfx = @(x,y) (-2*x)./(x.^2+rbf_c^2).^2*1./(y.^2+rbf_c^2);
% drbfy = @(x,y) (-2*y)./(y.^2+rbf_c^2).^2*1./(x.^2+rbf_c^2);

%rbf = @(t) exp(-rbf_c*t.^2); drbf = @(t) -2*rbf_c*exp(-rbf_c*t.^2);

bf = 0;
if bf == 1
    A_rbf = rbf((Xt_rbf-Xtau_rbf),(Yt_rbf-Ytau_rbf));
else
    A_rbf = rbf(R_rbf);
end
A_rbf = (A_rbf+A_rbf')/2;

rbf_mat = A_rbf;

%rbf_rhs = [rbf_f;zeros(p_basis_size,1)];
rbf_rhs = rbf_f';
%rbf_coeffs = rbf_mat\rbf_lhs;
rbf_coeffs = lsqminnorm(rbf_mat,rbf_rhs);
%poly_coefs = rbf_coeffs((end-p_basis_size+1):end);

points = mu_min:(mu_max-mu_min)/999:mu_max;

% rbf_fun = zeros(length(points),1);
% for i = 1:length(points)
%     point = points(i); 
%     for j = 1:length(rbf_points)
%         distr = norm(point-rbf_points(j)); 
%         rbf_fun(i) = rbf_fun(i)+rbf_coeffs(j)*rbf(distr);
%     end    
% end

rbf_fun = RBF(points,rbf,rbf_coeffs,rbf_points);
fun = @(x) RBF(x,rbf,rbf_coeffs,rbf_points);

lb = mu_min;             % Set a lower bound for the function.
ub = mu_max;          % Set an upper bound for the function.
nn = 1000;
x = NaN*ones(nn,1);             % Initializes x.
starting_points=linspace(lb,ub,nn);
for i=1:nn
        % Look for the zeros in the function's current window.
        x(i)=fzero(fun, starting_points(i));
end

mu_R = uniquetol(x);
%hold on; 
plot(points,rbf_fun,'color',"#7E2F8E",'LineWidth',3)
xline(mu_R,'color',"#000000","linewidth",2);
xlim([mu_min mu_max])
yline(0)
xlabel("$\mu$","Interpreter","latex")
ylabel("RBF","Interpreter","latex")

xticks(mu_R)
xtickangle(75)

title(join(["Robin-Laplace, $\sigma$ = ",sig],""),'Interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',20)

%sgtitle(join([curve_param,", RBF interpolating function and zeros"],""),'interpreter','latex','fontsize',20)

filename = join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/LapDNRev_minmu_',num2str(min(mus)),'maxmu_',num2str(max(mus)),curve_param,'_RBFzeros.png'],"");
filename = strrep(filename,"\","");
exportgraphics(gcf,filename,'Resolution',300)
%savefig(strrep(filename,".png",".fig"))

toc