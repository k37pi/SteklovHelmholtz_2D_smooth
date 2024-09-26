clear; close all; 
%% inputs -----------------------------------------------
a = 0; b = 2*pi; 
Mu = 2.5; % hopefully mu^2 is not a D/N ev 
N = 350;
%Ns = [50,100,150,200,250]; % N = half the number of intervals, set more than mu^2
%legtext = zeros(length(Ns),1);
Curve_name = ["coolpolar","coolpolar_1d"];%["ellipse","acorn"]; % out, in
if Mu > 0 
    mu = Mu;
else
    mu = -1i*Mu; % or +i Mu ?
end    
% % Epsilon test loop
% Disk_inr = 0.001:0.0001:0.01; %10.^(-(1:10));% 0.0001:0.0001:0.001;%
% eps_test = zeros(length(Disk_inr),2); 
% eps_test(:,1) = Disk_inr;
% %eps_test(:,2) = (2-mu^2)*pi./(2*pi*(1+Disk_inr));
% 
% for dinr = 1:length(Disk_inr)

% curve params
disk_R = [3,0.1]; % also ellipse a, diskdef_k, epps_egg
%disk_R = [1,Disk_inr(dinr)];
ell_b = [4,3];
acp = 0.55;
coolp_ab = [0.5000   -0.5000    0.1252    0.5000 0.2 0.37 0.45 0.6];
coolp_a0 = 1;

hom1 = 1.5;%3.44;
hom2 = 1;

len1 = 0;

% inside grid 
m = 250;
%--------------------------------------------------------
% for i = 1:length(Ns)
% N = Ns(i);
%% curve points and parametrization ---------------------
dt = (b-a)/(2*N); %cutoff = 0.1; 
t = a:dt:b; % 2N+1 points   

for dr = 1:2
disk_r = disk_R(dr); 
curve_name = Curve_name(dr);
switch curve_name
    case "circle"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,disk_r);
        curve_param = join([curve_name,disk_r],"");
    case "circle2"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,disk_r);
        curve_param = join([curve_name,disk_r],"");    
    case "ellipse"
        ellipse_a = disk_R(dr);
        ellipse_b = ell_b(dr);
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,ellipse_a,ellipse_b);
        curve_param = join([curve_name,num2str(ellipse_a),'\_',num2str(ellipse_b)],"");
    case "coolcurve"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,cool_ab);
        %curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3)),'\_',num2str(cool_ab(4)),'\_'],"");
        curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3))],"");
    case "disk_def"
        diskdef_k = disk_r;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,diskdef_k);  
        curve_param = join([curve_name,diskdef_k],"");
    case "rectangle"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,square_a,square_b);  
        curve_param = join([curve_name,num2str(square_a),'\_',num2str(square_b)],"");
      case "egg"
        epps_egg = disk_r;  
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,epps_egg);  
        curve_param = join([curve_name,num2str(epps_egg)],"");
    case 'coolpolar_1d'
        coolp_ab = acp;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,acp,1);
    case "coolpolar"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,coolp_ab,1);
    otherwise
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1);
        curve_param = curve_name;
end 

if dr == 1
x_out = hom1*x; dx_out = hom1*dx; d2x_out = hom1*d2x; nx_out = nx; len_out = hom1*len;
else
x_in = hom2*x; dx_in = hom2*dx; d2x_in = hom2*d2x; nx_in = -nx; len_in = hom2*len;
end    
end
T_len = len_out+len_in;
clear x;


xi = linspace(min(x_out(1,:)),max(x_out(1,:)),m);
yi = linspace(min(x_out(2,:)),max(x_out(2,:)),m);
[X,Y] = meshgrid(xi,yi);
idx = inpolygon(X(:),Y(:),x_out(1,:),x_out(2,:)) ;
points_inside1 = transpose([X(idx),Y(idx)]);

non_idx = inpolygon(points_inside1(1,:),points_inside1(2,:),x_in(1,:),x_in(2,:));
points_inside1(:,non_idx) = [];

fig=figure();
hold on;
% for i = 1:(2*N)
% plot([x_in(1,i) x_out(1,ceil(2*N*0.82))],[x_in(2,i) x_out(2,ceil(2*N*0.82))],'-')
% end
plot(x_out(1,:),x_out(2,:),'-*','LineWidth',1.4); hold on; 
%text(-disk_R(1),0,'$M_1$','Interpreter','latex','FontSize',20,'Color','blue')

plot(x_in(1,:),x_in(2,:),'-*','LineWidth',1.4); 
%text(-disk_R(2),0,'$M_2$','Interpreter','latex','FontSize',20,'Color','red')

%quiver(x_out(1,1),x_out(2,1),nx_out(1,1),nx_out(2,1),'LineWidth',3,'MaxHeadSize',1)
quiver(x_out(1,87),x_out(2,87),nx_out(1,87),nx_out(2,87),'color',"#EDB120",'LineWidth',3,'MaxHeadSize',1)
quiver(x_in(1,1),x_in(2,1),nx_in(1,1),nx_in(2,1),'LineWidth',3,'MaxHeadSize',1)
quiver(x_in(1,1),x_in(2,1),-nx_in(1,1),-nx_in(2,1),'LineWidth',3,'MaxHeadSize',1)

leg = legend(["$M_1$","$M_2$","$\nu_1$","$\nu_2$","$-\nu_2$"],"Interpreter","latex",'AutoUpdate','off');

plot(points_inside1(1,:),points_inside1(2,:),'o','color','#868686')
quiver(x_out(1,87),x_out(2,87),nx_out(1,87),nx_out(2,87),'color',"#EDB120",'LineWidth',3,'MaxHeadSize',1)
quiver(x_in(1,1),x_in(2,1),nx_in(1,1),nx_in(2,1),'Color',"#77AC30",'LineWidth',3,'MaxHeadSize',1)
quiver(x_in(1,1),x_in(2,1),-nx_in(1,1),-nx_in(2,1),'Color',"#7E2F8E",'LineWidth',3,'MaxHeadSize',1)

% for i = 1:(2*N)
% plot([x_in(1,ceil(2*N*0.82)) x_out(1,i)],[x_in(2,ceil(2*N*0.82)) x_out(2,i)],'-','LineWidth',1.4)
% end

plot(x_out(1,:),x_out(2,:),'*','color',"#0072BD",'LineWidth',3); hold on; 
%text(-disk_R(1),0,'$M_1$','Interpreter','latex','FontSize',20,'Color','blue')

plot(x_in(1,:),x_in(2,:),'*','color',"#D95319",'LineWidth',2); 
%text(-disk_R(2),0,'$M_2$','Interpreter','latex','FontSize',20,'Color','red')

%quiver(x_out(1,1),x_out(2,1),nx_out(1,1),nx_out(2,1),'color',"#EDB120",'LineWidth',3,'MaxHeadSize',1)
axis equal;
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20)
axis off; 
filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus.eps'];
%print(fig,'-depsc','-painters',filename);
% ax = gcf;
% exportgraphics(ax,strrep(filename,"eps","png"),'Resolution',300) 
% 
% ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% 
% hgexport(ax, filename,...
%     hgexport('factorystyle'), 'Format', 'eps');
% 
% legend(["Outer circle","Inner circle","$\nu_1$","$\nu_2$","$-\nu_2$"],"Interpreter","latex","FontSize",20)
% xticks([-disk_R(1) -disk_R(2) 0 disk_R(2) disk_R(1)])
% xticklabels({'$-R_1$','$-R_2$','0','$R_2$' '$R_1$'})
% yticks([-disk_R(1) -disk_R(2) 0 disk_R(2) disk_R(1)])
% yticklabels({'$-R_1$','$-R_2$','0','$R_2$' '$R_1$'})
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',20)
% filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus.png'];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
% saveas(gcf,filename);

%% true solution -------------------------
if sum(Curve_name == ["circle" "circle"])==2 %curve_name == "circle"
bessj1 = @(n) besselj(n,mu*disk_R(1)); bessj2 = @(n) besselj(n,mu*disk_R(2));
bessy1 = @(n) bessely(n,mu*disk_R(1)); bessy2 = @(n) bessely(n,mu*disk_R(2));
Dbessj1 = @(n) (besselj(n-1,mu*disk_R(1))-besselj(n+1,mu*disk_R(1)))/2; 
Dbessj2 = @(n) (besselj(n-1,mu*disk_R(2))-besselj(n+1,mu*disk_R(2)))/2;
Dbessy1 = @(n) (bessely(n-1,mu*disk_R(1))-bessely(n+1,mu*disk_R(1)))/2; 
Dbessy2 = @(n) (bessely(n-1,mu*disk_R(2))-bessely(n+1,mu*disk_R(2)))/2;
quadA = arrayfun(@(n) bessj1(n).*bessy2(n)-bessj2(n).*bessy1(n),0:N);   
quadB = mu*arrayfun(@(n) bessj1(n).*Dbessy2(n)-Dbessj1(n).*bessy2(n)-Dbessj2(n).*bessy1(n)+bessj2(n).*Dbessy1(n),0:N);
quadC = mu^2*arrayfun(@(n) -Dbessj1(n).*Dbessy2(n)+Dbessj2(n).*Dbessy1(n),0:N);
quadD = sqrt(quadB.^2-4*quadA.*quadC);

sig1 = (-quadB - quadD)./(2*quadA); sig1(isnan(sig1))=[]; sig1(isinf(sig1))=[]; 
Sig1 = sig1; sig1_mat = [Sig1',(0:length(Sig1)-1)'];
sigs1 = [sig1(1),repmat(sig1(2:end),1,2)]; sigs1 = sort(sigs1);

sig2 = (-quadB + quadD)./(2*quadA); sig2(isnan(sig2))=[]; sig2(isinf(sig2))=[]; 
Sig2 = sig2; sig2_mat = [Sig2',(0:length(Sig2)-1)'];
sigs2 = [sig2(1),repmat(sig2(2:end),1,2)]; sigs2 = sort(sigs2);

sig_mat = [sig1_mat;sig2_mat];

sig = [sigs1,sigs2]; sig(isinf(sig)) = []; sig(isnan(sig)) = []; Sig = sig';
[sig,ind_sig] = sort(sig); sig = sig';

k = 1; bessord = sig1_mat(sig1_mat(:,1)==sig1(k),2);
coeff_mat = [sig1(k)*bessj1(bessord) - mu*Dbessj1(bessord),sig1(k)*bessy1(bessord) - mu*Dbessy1(bessord); ...
    sig1(k)*bessj2(bessord) + mu*Dbessj2(bessord),sig1(k)*bessy2(bessord) + mu*Dbessy2(bessord)];
ann_a = 1;
ann_b = -ann_a*coeff_mat(1,1)/coeff_mat(1,2);
%ann_b = -ann_a*coeff_mat(2,1)/coeff_mat(2,2)

true_sol_out = (ann_a*bessj1(bessord)+ann_b*bessy1(bessord)).*exp(1i*bessord*t);
true_sol_in = (ann_a*bessj2(bessord)+ann_b*bessy2(bessord)).*exp(1i*bessord*t);
% what1 = Sig1.^2.*quadA+Sig1.*quadB+quadC;
% what2 = Sig2.^2.*quadA+Sig2.*quadB+quadC;
end
%% meshgrids ----------------------------------
[Tau,T] = meshgrid(t);

[Xtau_out,Xt_out] = meshgrid(x_out(1,:));
[Ytau_out,Yt_out] = meshgrid(x_out(2,:));
[dXtau_out,dXt_out] = meshgrid(dx_out(1,:));
[dYtau_out,dYt_out] = meshgrid(dx_out(2,:)); 
[d2Xtau_out,d2Xt_out] = meshgrid(d2x_out(1,:));
[d2Ytau_out,d2Yt_out] = meshgrid(d2x_out(2,:));

[Xtau_in,Xt_in] = meshgrid(x_in(1,:));
[Ytau_in,Yt_in] = meshgrid(x_in(2,:));
[dXtau_in,dXt_in] = meshgrid(dx_in(1,:));
[dYtau_in,dYt_in] = meshgrid(dx_in(2,:)); 
[d2Xtau_in,d2Xt_in] = meshgrid(d2x_in(1,:));
[d2Ytau_in,d2Yt_in] = meshgrid(d2x_in(2,:));

t_scaled = (0:2*N).*T; % the first column is 0*t, second is 1*t, so on.
cos_vecs = cos(t_scaled); % l is fixed in columns, t_j is fixed in rows.
sin_vecs = sin(t_scaled);
%% point on outer circle wrt to outer circle
[A1,B1] = layer_pots(mu,Tau,T,Xtau_out,Xt_out,Ytau_out,Yt_out,dXtau_out,dXt_out,dYtau_out,dYt_out,d2Xt_out,d2Yt_out);

%% point on inner circle wrt to inner circle
[A2,B2] = layer_pots(mu,Tau,T,Xtau_in,Xt_in,Ytau_in,Yt_in,dXtau_in,dXt_in,dYtau_in,dYt_in,d2Xt_in,d2Yt_in);

%% point on outer circle (required points) wrt to inner circle (x)
x = x_in; dXtau = dXtau_in; dYtau = dYtau_in; 
required_points = x_out; dXt = dXt_out; dYt = dYt_out; 

[Xtauin,Xin] = meshgrid(x(1,:),required_points(1,:)); % mesh x coord of req. points and bdry points
[Ytauin,Yin] = meshgrid(x(2,:),required_points(2,:)); % mesh y coord of req. points and bdry points

Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin; % tau corresponds to the boundary under consideration
Rin = sqrt(Xindiff.^2+Yindiff.^2); % rows correspond to required point 
dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);
dpnorm_ratio = dp2norm./dp1norm;

X = dYt.*Xindiff-dXt.*Yindiff; % outward normal for outer boundary dotted with distance between point on both boundaries
XR_ratio = X./Rin;
A2_tilde = (1i*mu/2)*dpnorm_ratio.*XR_ratio.*besselh(1,mu*Rin);
A2_tilde = pi/N*ftp(A2_tilde,N,cos_vecs,sin_vecs);
A2_tilde(:,end) = []; A2_tilde(end,:) = [];
% A2_tilde_tp = ftp(A2_tilde,N,cos_vecs,sin_vecs);
% [A2_tilde(2,:);A2_tilde_tp(2,:)]

B2_tilde = (1i/2)*besselh(0,mu*Rin).*dp2norm; % make this a tp and then sum and scale
B2_tilde = pi/N*ftp(B2_tilde,N,cos_vecs,sin_vecs);
B2_tilde(:,end) = []; B2_tilde(end,:) = [];

%B2_tilde_tp = ftp(B2_tilde,N,cos_vecs,sin_vecs);
%[B2_tilde(2,:);B2_tilde_tp(2,:)]

%% point on inner circle (required points) wrt to outer circle (x)
x = x_out; dXtau = dXtau_out; dYtau = dYtau_out; 
required_points = x_in; dXt = dXt_in; dYt = dYt_in; 

[Xtauin,Xin] = meshgrid(x(1,:),required_points(1,:)); % mesh x coord of req. points and bdry points
[Ytauin,Yin] = meshgrid(x(2,:),required_points(2,:)); % mesh y coord of req. points and bdry points

Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin; % tau corresponds to the boundary under consideration
Rin = sqrt(Xindiff.^2+Yindiff.^2); % rows correspond to required point 
dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);
dpnorm_ratio = dp2norm./dp1norm;

X = dYt.*Xindiff-dXt.*Yindiff; % outward normal for outer boundary dotted with distance between point on both boundaries
XR_ratio = X./Rin;
A1_tilde = (1i*mu/2)*dpnorm_ratio.*XR_ratio.*besselh(1,mu*Rin);
A1_tilde = pi/N*ftp(A1_tilde,N,cos_vecs,sin_vecs);
A1_tilde(:,end) = []; A1_tilde(end,:) = [];
% A1_tilde_tp = ftp(A1_tilde,N,cos_vecs,sin_vecs);
% [A1_tilde(2,:);A1_tilde_tp(2,:)]

B1_tilde = (1i/2)*besselh(0,mu*Rin).*dp2norm; % make this a tp and then sum and scale
B1_tilde = pi/N*ftp(B1_tilde,N,cos_vecs,sin_vecs);
B1_tilde(:,end) = []; B1_tilde(end,:) = [];
% B1_tilde_tp = ftp(B1_tilde,N,cos_vecs,sin_vecs);
% [B1_tilde(2,:);B1_tilde_tp(2,:)]

%% Matrix system ------------------------------
I = eye(2*N);
A = [A1+I,A2_tilde; -A1_tilde,-A2+I];
B = [B1,B2_tilde;B1_tilde,B2];
B_fms = [B1,zeros(2*N);zeros(2*N),B2];

%% Steklov eigenstuff ---------------------
[eden1,evs1] = eig(A,B); 
%[~,evs1,eden1] = eig(I+A,B);
evs1 = diag(evs1);  
%evs1 = real(evs1); 

[eden_in,evs_in1] = eig(A2+I,B2);
evs_in = diag(evs_in1);

[eden_out,evs_out1] = eig(A1+I,B1);
evs_out = diag(evs_out1);

if sum(imag(evs1)) < 1e-9
    evs = real(evs1);    
    [evs, ind] = sort(evs);%,"descend");
    eden = eden1(:,ind); % eigen densities
   
    [~,blah,ev_reps] = uniquetol(real(evs)); 
    ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals

    evs_in = real(evs_in); 
    [evs_in, ind_in] = sort(evs_in);%,"descend");
    eden_in = eden_in(:,ind_in); % eigen densities
    
    evs_out = real(evs_out); 
    [evs_out, ind_out] = sort(evs_out);%,"descend");
    eden_out = eden_out(:,ind_out); % eigen densities
else
    [evs, ind] = sort(evs1);%,"descend");
    eden = eden1(:,ind); % eigen densities
    [~,blah,ev_reps] = uniquetol(abs(evs)); 
    ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals
end    
evecs = B*eden/2;
evecs1 = B_fms*eden/2;

figure()
plot(evs,'--o','LineWidth',2); hold on; 
ylabel('$\sigma_k$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
set(gca,"FontSize",20)

if sum(Curve_name == ["circle" "circle"])==2 %curve_name == "circle"
%min_size = 2*N;
min_size = min([length(sig), length(evs)]);%, N]);
min_size = 400;
err = abs(sig(1:min_size)-evs(1:min_size)); rel_err = err./abs(sig(1:min_size));

fig = figure(3);
subplot(2,2,1:2); plot(1:min_size,evs(1:min_size),'--o','LineWidth',2); hold on; 
plot(1:min_size,sig(1:min_size),'--*','LineWidth',1.4); hold on;
plot(ev_reps(1),evs(ev_reps(1)),'o','MarkerFaceColor', 'b');
legend('Approx.','True','Single');
ylabel('$\sigma_k$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
set(gca,"FontSize",20)

subplot(2,2,3); semilogy(1:length(err),err,'--o','LineWidth',1.5); hold on; 
xlabel('$k$','Interpreter','latex'); ylabel('$AE_k$','Interpreter','latex')
set(gca,"FontSize",20)

subplot(2,2,4); 
%loglog(1:length(err),rel_err,'--o','LineWidth',1.5);
semilogy(1:length(err),rel_err,'--o','LineWidth',1.5);
xlabel('$k$','Interpreter','latex'); ylabel('$RE_k$','Interpreter','latex')
set(gca,"FontSize",20)

sgtitle(join(["Eigenvalues and errors on annulus with radii $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu,", $N$ = ",N],""),"interpreter",'latex','FontSize',20)

filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'_errs.png'],"")];
filename = ["/home/kshitij/Desktop/SFU/conferences/algoritmy 2024/talk/figs/annulus1.eps"];
%filename = [join(['august talk/annular/annulus_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'_errs.png'],"")];

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% hgexport(ax, strrep(filename,".png",".eps"),...
%     hgexport('factorystyle'), 'Format', 'eps');

% exportgraphics(ax,filename,'Resolution',300) 
% savefig(ax,strrep(filename,'.png','.fig'))
% print(fig,'-depsc','-vector',strrep(filename,'.png','.eps'));

%saveas(gcf,filename);

    % figure(3)
    % %subplot(2,3,i)
    % plot(evs,'LineWidth',1.5)
    % hold on; plot(sig,'LineWidth',1.5)
    % ylabel('$\sigma_k$','Interpreter','latex')
    % xlabel('$k$','Interpreter','latex')
    % set(gca,"FontSize",20)
    % legend('approx.','actual');
    % 
    %legtext(i) = join(["$N$ = ",N],"");
    % figure(4)
    % %subplot(2,3,i)
    % plot(sig,'--*','LineWidth',1.2); hold on;
    % ylabel('$\sigma_k$','Interpreter','latex')
    % xlabel('$k$','Interpreter','latex')
    % title(join(["True $\sigma_k$, $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu],""),'Interpreter','latex')
    % set(gca,"FontSize",20)
    %legend('approx.','actual');
    %legend([Ns]);
% end
% leg = legend(num2str(Ns'));
% title(leg,"$N$",'Interpreter','latex')
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
% filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_true_evs_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'.png'],"")];
% saveas(gcf,filename);
% ax = gcf;
% ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% hgexport(ax, strrep(filename,".png",".eps"),...
%     hgexport('factorystyle'), 'Format', 'eps');

end

% eps_test(dinr,2) = real(evs(1));
% close all
%end % eps test over

% figure();
% semilogx(eps_test(:,1),eps_test(:,2),'--*','LineWidth',2)
% hold on;
% %semilogx(eps_test(:,1),eps_test(:,2),'-->','LineWidth',2)
% xlabel("$\epsilon$",'Interpreter','latex')
% ylabel("$\sigma_1$",'Interpreter','latex')
% title(join(["Annulus with $R_1 = 1$, $R_2 = \epsilon$, $\mu$ = ",mu],""),'Interpreter','latex')
% set(gca,'FontSize',20)
% %filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_mu_',num2str(mu),'_N_',num2str(N),'_radius_eps_test.png'],"")];
% filename = [join(['august talk/annular/annulus_mu_',num2str(mu),'_N_',num2str(N),'_radius_eps_test_',num2str(Disk_inr(1)),'_',num2str(Disk_inr(end)),'.png'],"")];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
% set(gca, 'xdir', 'reverse')
%legend(["$\sigma_{1,\epsilon}$","$(2-\mu^2)\pi/|M_\epsilon|$"],'Interpreter','latex')
% exportgraphics(gcf,filename,'Resolution',300)
% savefig(gcf,strrep(filename,'.png','.fig'))
% writematrix(eps_test,strrep(filename,'.png','.csv'))

%% solution in Omega on M ---------------------------
k = 2;
if sum(Curve_name == ["circle" "circle"])==2 %curve_name == "circle"
bessord = sig_mat(abs(sig_mat(:,1)-evs(k)) < 1e-9,2);
coeff_mat = [evs(k)*bessj1(bessord) - mu*Dbessj1(bessord),evs(k)*bessy1(bessord) - mu*Dbessy1(bessord); ...
    evs(k)*bessj2(bessord) + mu*Dbessj2(bessord),evs(k)*bessy2(bessord) + mu*Dbessy2(bessord)];
ann_a = 1;
ann_b = -ann_a*coeff_mat(1,1)/coeff_mat(1,2);
%ann_b = -ann_a*coeff_mat(2,1)/coeff_mat(2,2)

true_sol_out = (ann_a*bessj1(bessord)+ann_b*bessy1(bessord)).*exp(1i*bessord*t(1:end-1));
true_sol_in = (ann_a*bessj2(bessord)+ann_b*bessy2(bessord)).*exp(1i*bessord*t(1:end-1));

subspace(real(evecs(1:2*N,k)),[real(true_sol_out);imag(true_sol_out)]')
real_out_ang = subspace(real(evecs(1:2*N,k)),[real(true_sol_out);imag(true_sol_out)]')*180/pi;

subspace(imag(evecs(1:2*N,k)),[real(true_sol_out);imag(true_sol_out)]')
imag_out_ang = subspace(imag(evecs(1:2*N,k)),[real(true_sol_out);imag(true_sol_out)]')*180/pi;

subspace(real(evecs((2*N+1):end,k)),[real(true_sol_in);imag(true_sol_in)]')
real_in_ang = subspace(real(evecs((2*N+1):end,k)),[real(true_sol_in);imag(true_sol_in)]')*180/pi;

subspace(imag(evecs((2*N+1):end,k)),[real(true_sol_in);imag(true_sol_in)]')
imag_in_ang = subspace(imag(evecs((2*N+1):end,k)),[real(true_sol_in);imag(true_sol_in)]')*180/pi;
end

dp2norm_out = sqrt(dXtau_out.^2+dYtau_out.^2);
dp2norm_in = sqrt(dXtau_in.^2+dYtau_in.^2);
eden_fix = eden(:,k); eden_fix_out = eden_fix(1:2*N); eden_fix_in = eden_fix((2*N+1):end); 

for ins = 1:2
    if ins == 1
        [Xtauin,Xin] = meshgrid(x_out(1,1:end-1),points_inside1(1,:)); % mesh x coord of in points and bdry points
        [Ytauin,Yin] = meshgrid(x_out(2,1:end-1),points_inside1(2,:)); % mesh y coord of in points and bdry points
        Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
        Rin = sqrt(Xindiff.^2+Yindiff.^2);
        HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
        Recon_integrand = HankRin.*transpose(eden_fix_out).*dp2norm_out(k,(1:end-1)); % make this a tp and then sum and scale
        recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
        recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
        recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
        recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
        recint_sin = recint_sin1*sin_vecs(:,2:N)';
        recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
        %solution_fixed_point = 1i*pi/(4*N)*sum(Recon_integrand,2);
        solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
        %solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 
    else
        [Xtauin,Xin] = meshgrid(x_in(1,1:end-1),points_inside1(1,:)); % mesh x coord of in points and bdry points
        [Ytauin,Yin] = meshgrid(x_in(2,1:end-1),points_inside1(2,:)); % mesh y coord of in points and bdry points
        Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
        Rin = sqrt(Xindiff.^2+Yindiff.^2);
        HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
        Recon_integrand = HankRin.*transpose(eden_fix_in).*dp2norm_in(k,(1:end-1)); % make this a tp and then sum and scale
        recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
        recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
        recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
        recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
        recint_sin = recint_sin1*sin_vecs(:,2:N)';
        recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
        %solution_fixed_point = 1i*pi/(4*N)*sum(Recon_integrand,2);
        solution_fixed_point = solution_fixed_point+1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
        %solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1));       
    end    
end

% xlin = linspace(min(points_inside1(1,:)),max(points_inside1(1,:)),round(1*m));
% ylin = linspace(min(points_inside1(2,:)),max(points_inside1(2,:)),round(1*m));
xlin = linspace(min([x_out(1,:) x_in(1,:)]),max([x_out(1,:) x_in(1,:)]),round(1*m));
ylin = linspace(min([x_out(2,:) x_in(2,:)]),max([x_out(2,:) x_in(2,:)]),round(1*m));
[Xin2,Yin2] = meshgrid(xlin,ylin);

all_points = [points_inside1,x_out(:,1:end-1),x_in(:,1:end-1)];
all_sol = [solution_fixed_point;evecs(:,k)];
    
%Uin = griddata(points_inside1(1,:),points_inside1(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
Uin = griddata(all_points(1,:),all_points(2,:),all_sol,Xin2,Yin2,'cubic');

RealUin = real(Uin);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
%RealUin = RealUin/max(abs(solution_fixed_point));
ImagUin = imag(Uin);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');
%ImagUin = ImagUin/max(abs(solution_fixed_point));

%%{
bdrypts = inpolygon(Xin2,Yin2,x_in(1,:),x_in(2,:));
Xin2(bdrypts) = NaN; Yin2(bdrypts) = NaN;
bdrypts2 = inpolygon(Xin2,Yin2,x_out(1,:),x_out(2,:));
Xin2(~bdrypts2) = NaN; Yin2(~bdrypts2) = NaN;
%}

if sum(Curve_name == ["circle" "circle"])==2 %curve_name == "circle"
Rin2 = abs(mu)*sqrt(Xin2.^2+Yin2.^2); 
BessRin2 = ann_a*besselj(bessord,Rin2)+ann_b*bessely(bessord,Rin2);
cos_in = cos(bessord*atan2(Yin2,Xin2));
sin_in = sin(bessord*atan2(Yin2,Xin2));
Z_real_ac = BessRin2.*cos_in;
Z_imag_ac = BessRin2.*sin_in;
Z1 = Z_real_ac+1i*Z_imag_ac; Z2 = Z_real_ac-1i*Z_imag_ac; 

%%{
app_sol_vec = reshape(Uin,[],1); 
ac_sol1 = reshape(Z1,[],1); ac_sol1(isnan(app_sol_vec)) = [];
app_sol_vec(isnan(app_sol_vec)) = [];
app_sol_vec(isnan(ac_sol1)) = [];
ac_sol1(isnan(ac_sol1)) = [];
%}

real_inside_ang = subspace(real(app_sol_vec), [real(ac_sol1) imag(ac_sol1)])*180/pi;
imag_inside_ang = subspace(imag(app_sol_vec), [real(ac_sol1) imag(ac_sol1)])*180/pi;

figure()
subplot(2,2,1)
surf(Xin2,Yin2,RealUin);
shading interp; view(0,90); colormap(magma); colorbar
title('$Real(u^N_k)$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Real$(u)$','Interpreter','latex');
iptsetpref('ImshowBorder','tight');
ylim([-disk_R(1) disk_R(1)])
xlim([-disk_R(1) disk_R(1)])
axis equal

subplot(2,2,2)  
surf(Xin2,Yin2,real(Z1));
shading interp; view(0,90); colormap(magma); colorbar
title('$Real(u_k)$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Real$(u)$','Interpreter','latex');
iptsetpref('ImshowBorder','tight');
ylim([-disk_R(1) disk_R(1)])
xlim([-disk_R(1) disk_R(1)])
axis equal

subplot(2,2,3)
surf(Xin2,Yin2,ImagUin);
shading interp; view(0,90); colormap(magma); colorbar
title('$Imag(u^N_k)$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Real$(u)$','Interpreter','latex');
ylim([-disk_R(1) disk_R(1)])
xlim([-disk_R(1) disk_R(1)])
axis equal

subplot(2,2,4)  
surf(Xin2,Yin2,imag(Z1));
shading interp; view(0,90); colormap(magma); colorbar
title('$Imag(u_k)$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Real$(u)$','Interpreter','latex');
iptsetpref('ImshowBorder','tight');
ylim([-disk_R(1) disk_R(1)])
xlim([-disk_R(1) disk_R(1)])
axis equal

sgtitle({ join(["Annulus, $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu,", $N$ = ",N,", $k$ = ",k],"")...
    join(["Angle between $u_k$ and $Real(u^N_k)$ = ",real_inside_ang],"")...
    join(["Angle between $u_k$ and $Imag(u^N_k)$ = ",imag_inside_ang],"")...
},'Interpreter','latex','FontSize',18)
set(findobj(gcf,'type','axes'),'FontSize',18)
filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_evecs_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'_k',num2str(k),'_inside.png'],"")];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
exportgraphics(gcf,filename,'Resolution',300)  
savefig(gcf,strrep(filename,".png",".fig"))
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% hgexport(ax, strrep(filename,".png",".eps"),...
%     hgexport('factorystyle'), 'Format', 'eps');


figure()
subplot(2,2,1)
plot3(x_out(1,1:end-1),x_out(2,1:end-1),real(evecs(1:2*N,k)))
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),real(evecs((2*N+1):end,k)))

subplot(2,2,3)
plot3(x_out(1,1:end-1),x_out(2,1:end-1),imag(evecs(1:2*N,k)))
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),imag(evecs((2*N+1):end,k)))

subplot(2,2,2)
plot3(x_out(1,1:end-1),x_out(2,1:end-1),real(true_sol_out))
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),real(true_sol_in))

subplot(2,2,4)
plot3(x_out(1,1:end-1),x_out(2,1:end-1),imag(true_sol_out))
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),imag(true_sol_in))


figure()
subplot(2,2,1)
plot(t(1,1:end-1),real(evecs(1:2*N,k)),'LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Real(u^N_k(R_1 t))$",'Interpreter','latex')

subplot(2,2,2)
plot(t(1:end-1),real(true_sol_out),'r','LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Real(u_k(R_1 t))$",'Interpreter','latex')

subplot(2,2,3)
plot(t(1,1:end-1),imag(evecs(1:2*N,k)),'LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Imag(u^N_k(R_1 t))$",'Interpreter','latex')

subplot(2,2,4)
plot(t(1:end-1),imag(true_sol_out),'r','LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Imag(u_k(R_1 t))$",'Interpreter','latex')

sgtitle({ join(["Annulus, $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu,", $N$ = ",N,", $k$ = ",k],"")...
    join(["Angle between $u_k$ and $Real(u^N_k(R_1 t))$ = ",real_out_ang],"")...
    join(["Angle between $u_k$ and $Imag(u^N_k(R_1 t))$ = ",imag_out_ang],"")...
},'Interpreter','latex','FontSize',23)
set(findobj(gcf,'type','axes'),'FontSize',23)
filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_evecs_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'_k',num2str(k),'_onM1.png'],"")];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
exportgraphics(gcf,filename,'Resolution',300)  
savefig(gcf,strrep(filename,".png",".fig"))
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% hgexport(ax, strrep(filename,".png",".eps"),...
%     hgexport('factorystyle'), 'Format', 'eps');

figure()
subplot(2,2,1)
plot(t(1:end-1),real(evecs((2*N+1):end,k)),'LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Real(u^N_k(R_2 t))$",'Interpreter','latex')

subplot(2,2,2)
plot(t(1:end-1),real(true_sol_in),'r','LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Real(u_k(R_2 t))$",'Interpreter','latex')

subplot(2,2,3)
plot(t(1:end-1),imag(evecs((2*N+1):end,k)),'LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Imag(u^N_k(R_2 t))$",'Interpreter','latex')

subplot(2,2,4)
plot(t(1:end-1),imag(true_sol_in),'r','LineWidth',2)
ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Imag(u_k(R_2 t))$",'Interpreter','latex')

sgtitle({ join(["Annulus, $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu,", $N$ = ",N,", $k$ = ",k],"")...
    join(["Angle between $u_k$ and $Real(u^N_k(R_2 t))$ = ",real_in_ang],"")...
    join(["Angle between $u_k$ and $Imag(u^N_k(R_2 t))$ = ",imag_in_ang],"")...
},'Interpreter','latex','FontSize',23)

set(findobj(gcf,'type','axes'),'FontSize',23)
filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_evecs_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'_k',num2str(k),'_onM2.png'],"")];

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
exportgraphics(gcf,filename,'Resolution',300)  
savefig(gcf,strrep(filename,".png",".fig"))
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
% hgexport(ax, strrep(filename,".png",".eps"),...
%     hgexport('factorystyle'), 'Format', 'eps');

else
figure()
subplot(1,1,1)
surf(Xin2,Yin2,RealUin); 
%surf(Xin2,Yin2,abs(Uin)); 
shading interp; 
view(0,90); colormap(magma); colorbar
hold on;
plot3(x_out(1,1:end-1),x_out(2,1:end-1),real(evecs(1:2*N,k)),'LineWidth',2)
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),real(evecs((2*N+1):end,k)),'LineWidth',2)
% plot3(x_out(1,1:end-1),x_out(2,1:end-1),abs(evecs(1:2*N,k)),'LineWidth',2)
% hold on
% plot3(x_in(1,1:end-1),x_in(2,1:end-1),abs(evecs((2*N+1):end,k)),'LineWidth',2)
title('$Real(u^N_k)$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
% das = daspect; das(1) = das(2);
% set(gca,'DataAspectRatio',das);
axis equal

% subplot(1,2,2)
% surf(Xin2,Yin2,ImagUin); 
% %surf(Xin2,Yin2,abs(Uin)); 
% shading interp; 
% view(0,90); colormap(magma); colorbar
% hold on;
% plot3(x_out(1,1:end-1),x_out(2,1:end-1),imag(evecs(1:2*N,k)),'LineWidth',2)
% hold on;
% plot3(x_in(1,1:end-1),x_in(2,1:end-1),imag(evecs((2*N+1):end,k)),'LineWidth',2)
% % plot3(x_out(1,1:end-1),x_out(2,1:end-1),abs(evecs(1:2*N,k)),'LineWidth',2)
% % hold on;
% % plot3(x_in(1,1:end-1),x_in(2,1:end-1),abs(evecs((2*N+1):end,k)),'LineWidth',2)
% title('$Real(u^N_k)$','Interpreter','latex')
% xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
% % das = daspect; das(1) = das(2);
% % set(gca,'DataAspectRatio',das);
% axis equal


figure()
subplot(1,1,1)
surf(Xin2,Yin2,abs(Uin)); 
%surf(Xin2,Yin2,abs(Uin)); 
shading interp; 
view(0,90); colormap(magma); colorbar
hold on;
plot3(x_out(1,1:end-1),x_out(2,1:end-1),abs(evecs(1:2*N,k)),'LineWidth',2)
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),abs(evecs((2*N+1):end,k)),'LineWidth',2)
% plot3(x_out(1,1:end-1),x_out(2,1:end-1),abs(evecs(1:2*N,k)),'LineWidth',2)
% hold on;
% plot3(x_in(1,1:end-1),x_in(2,1:end-1),abs(evecs((2*N+1):end,k)),'LineWidth',2)
title('$|u^N_k|$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
% das = daspect; das(1) = das(2);
% set(gca,'DataAspectRatio',das);
axis equal

figure()
subplot(2,1,1)
plot3(x_out(1,1:end-1),x_out(2,1:end-1),real(evecs(1:2*N,k)))
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),real(evecs((2*N+1):end,k)))

subplot(2,1,2)
plot3(x_out(1,1:end-1),x_out(2,1:end-1),imag(evecs(1:2*N,k)))
hold on;
plot3(x_in(1,1:end-1),x_in(2,1:end-1),imag(evecs((2*N+1):end,k)))


figure()
subplot(2,1,1)
plot(t(1,1:end-1),real(evecs(1:2*N,k)))
%ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Real(u^N_k(R_1 t))$",'Interpreter','latex')

subplot(2,1,2)
plot(t(1,1:end-1),imag(evecs(1:2*N,k)))
%ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Imag(u^N_k(R_1 t))$",'Interpreter','latex')

sgtitle(join(["Annulus, $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu,", $N$ = ",N,", $k$ = ",k],""),'Interpreter','latex','FontSize',18)
set(findobj(gcf,'type','axes'),'FontSize',18)
filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_evecs_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'out_k',num2str(k),'.png'],"")];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
%saveas(gcf,filename);

figure()
subplot(2,1,1)
plot(t(1:end-1),real(evecs((2*N+1):end,k)))
%ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Real(u^N_k(R_2 t))$",'Interpreter','latex')

subplot(2,1,2)
plot(t(1:end-1),imag(evecs((2*N+1):end,k)))
%ylim([-1 1])
xlabel("$t$",'Interpreter','latex')
ylabel("$Imag(u^N_k(R_2 t))$",'Interpreter','latex')

sgtitle(join(["Annulus, $R_1$ = ",disk_R(1),", $R_2$ = ",disk_R(2),", $\mu$ = ",mu,", $N$ = ",N,", $k$ = ",k],"")...
,'Interpreter','latex','FontSize',18)

set(findobj(gcf,'type','axes'),'FontSize',18)
filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/annulus_evecs_mu_',num2str(mu),'_N_',num2str(N),'_',num2str(disk_R(1)),'_',num2str(disk_R(2)),'in_k',num2str(k),'.png'],"")];

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
%saveas(gcf,filename);
 
end
    
evs(1) < (2-mu^2)*pi/(T_len)
evs(1) < (2-mu^2)*pi/2*(1/(len_out)+1/(len_in))

figure()
plot(evs(1:N/2),'--*','LineWidth',2)
hold on
plot(evs_in(1:N/2),'--*','LineWidth',2)
plot(evs_out(1:N/2),'--*','LineWidth',2)
legend(["Annuluar","In","Out"])