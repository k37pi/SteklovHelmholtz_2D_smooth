tic
clear; close all; 
%% inputs -----------------------------------------------
a = 0; b = 2*pi;  
Mu = 5;%.42275;2.41355;1.96702;1.1698; 1.164; %1.693;%3.831;%8.77;%.5678; % hopefully mu^2 is not a D/N ev 
N = 350; % N = half the number of intervals, set more than mu^2
curves_list = ["" + ...
    "Circle"; "Ellipse"; "kite";"coolcurve";"disk_def";
    "acorn";"star";"rectangle";"egg";"opt_disk";"butterfly";"cavity";"bunimovich"];%;"peanut"];
% cool curve = (a1 cos (b1 t) + a2 sin (b2 t))ang2a[cost, sint] 

% circle = 1, ellipse = 2,  kite = 3, coolcurve = 4, 
% disk_def = 5, acorn = 6 , star = 7, square = 8, egg = 9
% opt_disk = 10, butterfly = 11, cavity = 12, bunimovich = 13
curve_name = curves_list(11); curve_name = lower(curve_name);
rho_wt = @(t) 1;%sin(3*t).^2+1;%exp(cos(t));

if Mu > 0 
    mu = Mu;
else
    mu = -1i*Mu; % or +i Mu ?
end    

% curve params
hom = 1; % domain scaling
disk_r = 7;%/(2*pi);%3/(2*pi);%0.208;%3;
ellipse_a = 21; ellipse_b = 14;
%cool_ab = [2,2,3];
square_a = 1; square_b = 2;
epps_egg = 0.1;
aa = -1/4; bb = 1/4;
cool_ab =  [0.2500   -0.0388    0.1634   -0.0473];%aa+(bb-aa)*rand(4,1);%randi([1,70],2,1); 
%cool_ab = [cool_ab;randi([1,100],1,1)]; 
% cool_ab = [cool_ab;randi([cool_ab(1),100],1,1)]; 
% c1 = (cool_ab(3)^2-cool_ab(1)^2)/(4*pi^2*cool_ab(1)^2);
% b1 = sqrt((1+sqrt(1+4*c1))/2);
% cool_ab(2) = randi([floor(0.5*(cool_ab(1)+cool_ab(3))),100],1,1);
%cool_ab = [33;3;68];
diskdef_k = 1.5;%.4448;%-192.1716;%-5.30538947734613e-06; 
len1 = 0;

% inside grid
m = 300;
%--------------------------------------------------------

%% curve points and parametrization ---------------------
dt = (b-a)/(2*N); %cutoff = 0.1; 
t = a:dt:b; % 2N+1 points
% t1 = arrayfun(@(n) pi*(1+cos((2*n)*pi/(4*N))), 0:2*N);
% t=sort(t1);
%t1 = [t1 0]; t1 = sort(t1); t1(end) = 2*pi; t = t1;

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
        while(~isempty(selfintersect(x(1,1:end-1),x(2,1:end-1))))
            cool_ab = aa+(bb-aa)*rand(4,1);
            [x,dx,d2x,nx,len] = curve(curve_name,t,len1,cool_ab);
            curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3))],"");
        end    

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

if curve_name ~= "circle" && len1 ~= 1
    x = hom*x; dx = hom*dx; d2x = hom*d2x; len = hom*len;
    nx = [dx(2,:);-dx(1,:)];
    nx = nx./sqrt(dx(1,:).^2+dx(2,:).^2);
end

max(dx(1,:).^2+dx(2,:).^2)/min(dx(1,:).^2+dx(2,:).^2)

dom_area = polyarea(x(1,:)',x(2,:)');
kappa = (dx(1,1:end-1).*d2x(2,1:end-1)-d2x(1,1:end-1).*dx(2,1:end-1))./(dx(1,1:end-1).^2+dx(2,1:end-1).^2).^1.5;
kappa = sum(kappa)*pi/N; kappa = kappa/(2*pi);

xi = linspace(min(x(1,:)),max(x(1,:)),m);
yi = linspace(min(x(2,:)),max(x(2,:)),m);
[X,Y] = meshgrid(xi,yi) ;
idx = inpolygon(X(:),Y(:),x(1,:),x(2,:)) ;
points_inside = transpose([X(idx),Y(idx)]);

points_inside_close = 0.99*x(:,1:end-1);

points_inside = [points_inside,points_inside_close];

[Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
[Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

fig=figure(); 
hold on; 
plot(x(1,:),x(2,:),'-',"LineWidth",2.5); 
quiver(x(1,1),x(2,1),nx(1,1),nx(2,1),'LineWidth',2.5,'MaxHeadSize',0.5,'color',"#77AC30")
% for i = 1:(2*N)
% plot(x(1,[ceil(2*N*0.82) i]),x(2,[ceil(2*N*0.82) i]),'-')
% end
% plot(x(1,:),x(2,:),'--*',"color","#0072BD","LineWidth",2.4); 
% plot(points_inside(1,:),points_inside(2,:),'o');
text(-0.9,x(2,N),'$M$','Interpreter','latex','FontSize',30,'Color',	"#0072BD")
%text(0.5,0.5,'$\Omega$, $n_{RI}$','Interpreter','latex','FontSize',30)
text(0.5,0.5,'$\Omega$','Interpreter','latex','FontSize',30)
hold on;
quiver(x(1,1),x(2,1),nx(1,1),nx(2,1),'LineWidth',2.5,'MaxHeadSize',0.5,'color',"#77AC30")
text(3,0.14,'$\nu$','Interpreter','latex','FontSize',30,'color',"#77AC30")
% plot(x(1,3*N/2),x(2,3*N/2),'o','LineWidth',5,'color',"#A2142F")
% text(-1.33,x(2,3*N/2),"$(z_1(t_i),z_2(t_i))$",'FontSize',35,'Interpreter','latex','color',"#A2142F",'FontWeight','bold')
axis equal;
axis off
leg = legend(["$M = \partial \Omega$","$\nu$ = outer unit normal"],"Interpreter","latex",'AutoUpdate','off');
%plot(x(1,ceil(2*N*0.82)),x(2,ceil(2*N*0.82)),'*',"color","#000000",'LineWidth',12)
set(leg,'FontSize',20);

filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/domain.eps'];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
ax = gcf;
% exportgraphics(ax,filename,'Resolution',300) 
% saveas(gcf,filename);
print(fig,'-depsc','-painters',filename);

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
%toc
% ----------------------------------------
%% L, M ---------------------------------------
tic
dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);
dpnorm_ratio = dp2norm./dp1norm;

ind_ratios = R~=0;

%mu/2*max(XR_ratio(ind_ratios))*max(dpnorm_ratio(ind_ratios))

logsin_term = log(4*sin((T-Tau)/2).^2);
cL = (1i*mu/2)*dpnorm_ratio;
cL1 = -mu/(2*pi)*dpnorm_ratio;
L = cL.*XR_ratio.*besselh(1,mu*R); % note that besselh(1,0) blows up so when R = 0 it's useless.
L(end,1) = 0; L(1,end) = 0;

L1 = cL1.*XR_ratio.*besselj(1,mu*R);
% if Mu < 0
%     L1 = 1i*L1;
% end    

L2 = L-L1.*logsin_term;
L2(end,1) = 0; L2(1,end) = 0;
% if curve_name == "rectangle"
%     L(isnan(L)) = 0; 
%     L2(isnan(L2)) = 0; 
% else
    L(isnan(L)) = diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2)); 
    L2(isnan(L2)) = diag(L); 
% end
    
L1(isnan(L1)) = 0;


M = (1i/2)*dp2norm.*besselh(0,mu*R);
M(1,end) = 0; M(end,1) = 0;
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
M2(1,end) = 0; M2(end,1) = 0;
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
tic

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
    % %%toc
    % inds = (500*(ncount-1)+1):(N-1);
    % rep_diff1 = rep_diff1(:,:,1:length(inds)); % create 3d matrix with copies of diff
    % % rep_diff = rep_diff.*reshape(inds,1,1,[]); % multiply each matrix with next integer
    % % rep_diff_cos = cos(rep_diff); 
    % % rep_diff_cos_scaled = rep_diff_cos./reshape(inds,1,1,[]);
    % rep_diff_cos_scaled = cos(rep_diff1.*reshape(inds,1,1,[]))./reshape(inds,1,1,[]);
    % RN_mat = RN_mat -2*pi*sum(rep_diff_cos_scaled,3)/N - pi*cos(N*(diff_t_tau))/N^2;


    RN_mat = cos(diff_t_tau);
    for k1 = 2:(N-1)
        RN_mat = RN_mat+cos(k1*diff_t_tau)/k1; 
    end    
    RN_mat = -2*pi*RN_mat/N - pi*cos(N*(diff_t_tau))/N^2;
% end
%toc

%-------------------------------------------

%% Matrices ------------------------------
A = RN_mat.*L1_tp+pi/N*L2_tp; % (A + I) is the effect of K'
A(:,end) = []; A(end,:) = [];

I = diag(ones(2*N,1));

B = RN_mat.*M1_tp+pi/N*M2_tp; % B is the effect of the S 
B(:,end) = []; B(end,:) = [];
%B(1:2*N,1:2*N) = 0;

Rho_wt = rho_wt(t(1:end-1))';
B_wt = Rho_wt.*B;
% [A1,B1] = layer_pots(mu,Tau,T,Xtau,Xt,Ytau,Yt,dXtau,dXt,dYtau,dYt,d2Xt,d2Yt);
% sum(sum(A1==A)) == 4*N^2
% sum(sum(B1==B)) == 4*N^2

%% Steklov eigenstuff --------------------------------
%evs = eigs(I+A,B,5,'smallestabs');
[eden1,evs1] = eig(I+A,B_wt); 
%[eden1,evs1] = eig((B_wt)\(I+A)); 
evs1 = diag(evs1);
%evs1 = real(evs1); 

if sum(imag(evs1)) < 1e-9
    evs = real(evs1);    
    [evs, ind] = sort(evs);%,"descend");
    eden = eden1(:,ind); % eigen densities
    [~,blah,ev_reps] = uniquetol(real(evs)); 
    ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals
else
    [evs, ind] = sort(evs1);%,"descend");
    eden = eden1(:,ind); % eigen densities
    [~,blah,ev_reps] = uniquetol(abs(evs)); 
    ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals
end    
evecs = B_wt*eden/2; %disk_r
normder_evecs = (A+I)*eden/2;
% only for circle, be careful
if Mu > 0
    evs_ac = arrayfun(@(n) (mu)*0.5*(besselj(n-1,(mu*len/(2*pi)))-besselj(n+1,(mu*len/(2*pi))))./besselj(n,(mu*len/(2*pi))),0:N);
    %evs_ac = arrayfun(@(n) (mu)*0.5*(besselj(n-1,(mu*disk_r))-besselj(n+1,(mu*disk_r)))./besselj(n,(mu*disk_r)),0:N);
    % Big_besselj = @(n,mu) 1/sqrt(2*pi*n)*(exp(1)*mu/(2*n))^n;
    % evs_ac_big = arrayfun(@(n) (mu)*0.5*(Big_besselj(n-1,(mu))-Big_besselj(n+1,(mu)))./Big_besselj(n,(mu)),0:N);
    % evs_ac_big = [evs_ac_big(1),repmat(evs_ac_big(2:end),1,2)]; evs_ac_big = sort(evs_ac_big);
else
    evs_ac = arrayfun(@(n) abs(Mu)*0.5*(besseli(n-1,abs(Mu)*len/(2*pi))+besseli(n+1,abs(Mu)*len/(2*pi)))./besseli(n,abs(Mu)*len/(2*pi)),0:N);
end    
evsac_mat = [evs_ac',(0:N)'];
evs_ac = [evs_ac(1),repmat(evs_ac(2:end),1,2)]; [evs_ac,ind2] = sort(evs_ac);
evs_ac = evs_ac';

if curve_name == "circle"
    figure()
    subplot(2,2,1:2); plot(1:N,evs(1:N),'--o','LineWidth',2); hold on; 
    plot(1:N,evs_ac(1:N),'--*','LineWidth',2); hold on;
    plot(ev_reps(1),evs(ev_reps(1)),'o','MarkerFaceColor', 'b');
    % if N > 70  
    %     plot(evs_ac_big,'diamond');
    %     legend('approx.','actual','single','actual asymp','Location','north')
    % else
        legend('Approx.','True','Single');%,'Location','north')
    % end    
    %title(['$\mu=$ ',num2str(mu),', \# of points = ',num2str(2*N+1)],'Interpreter','latex')
    ylabel('$\sigma_k$','Interpreter','latex')
    xlabel('$k$','Interpreter','latex')
    %filename = ['figs/',curve_name,'_mu',num2str(mu),'_N',num2str(N),'.png'];
    %filename = join(filename,'');
    %saveas(gcf,filename);
    
    % figure()
    % plot(real(evecs(:,ev_reps(1))),'--o'); 
    % legend('approx.','actual','Location','north')
    % title(['$\mu=$ ',num2str(mu),', \# of points = ',num2str(2*N+1)],'Interpreter','latex')
    % ylabel('$\sigma$','Interpreter','latex')
    % xlabel('N')
    % filename = ['figs/',curve_name,'_mu',num2str(mu),'_N',num2str(N),'_evec',num2str(ev_reps(1)),'.png'];
    % saveas(gcf,filename);
    % evs(ev_reps)
    
    % plot_ind = 2; %ev_reps(1);
    % figure()
    % plot3(real(evecs(:,plot_ind)),imag(evecs(:,plot_ind)),t(1:end-1));
    % xlabel('Re(u(t))')
    % ylabel('Im(u(t))')
    % zlabel('t')

%% Errors (only for circle) ----------------------------
evs_ac_inf = find(evs_ac==Inf);
if ~isempty(evs_ac_inf)
    evs_ac_inf1 = evs_ac_inf(1);
    err = abs(evs(1:(evs_ac_inf1-1))-evs_ac(1:(evs_ac_inf1-1)));
    rel_err = abs(err./evs_ac(1:(evs_ac_inf1-1)));
else
    err = abs(evs-evs_ac(1:length(evs)));
    rel_err = err./abs(evs_ac(1:length(evs)));
end
%figure()
min_size = min([N,length(err)]);
%subplot(1,2,1); plot(err(1:end-round(N/5)),'--o'); hold on; 
subplot(2,2,3); plot(1:min_size,err(1:min_size),'--o','LineWidth',2); hold on; 
%title('abs. err.')
xlabel('$k$','Interpreter','latex'); ylabel('$AE_k$','Interpreter','latex')
%subplot(1,2,2); loglog(rel_err(1:end-round(N/10)),'--o');
subplot(2,2,4); semilogy(1:min_size,rel_err(1:min_size),'--o','LineWidth',2);
xlabel('$k$','Interpreter','latex'); ylabel('$RE_k$','Interpreter','latex')
%title('rel. err.')
%sgtitle(['$\mu=$ ',num2str(mu),', \# of points = ',num2str(2*N+1),' minus end'],'Interpreter','latex')
sgtitle(['Eigenvalues and errors on disk with $R=$ ',num2str(len/(2*pi)),', $\mu=$ ',num2str(mu),', $N=$ ',num2str(N)],'Interpreter','latex','fontsize',20) %,', \# of points = ',num2str(2*N+1)
fig_curvename = strrep(curve_param, '.', 'dot');
%filename = ['figs/',fig_curvename,'_mu',num2str(mu),'_N',num2str(N),'_errs.svg'];
filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/',fig_curvename,'_mu',num2str(mu),'_N',num2str(N),'_errs.png'];
filename = join(filename,'');
set(findobj(gcf,'type','axes'),'FontSize',20)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
exportgraphics(gcf,filename,'Resolution',300) 
saveas(gcf,filename);
else
    fig = figure();
    if sum(imag(evs)) == 0 
        plot(evs,'--o'); hold on; 
        plot(ev_reps,evs(ev_reps),'o','MarkerFaceColor', 'b');
        legend('approx.','single','Location','north')
        title(join([curve_param,', $\mu=$ ',num2str(mu),', \# of points = ',num2str(2*N+1)],""),'Interpreter','latex')
        ylabel('$\sigma$','Interpreter','latex')
        xlabel('N')
    else
        subplot(1,2,1); plot(real(evs),'--o'); hold on; 
        ylabel("Real($\sigma$)",'Interpreter','latex')
        subplot(1,2,2); plot(imag(evs),'--o'); hold on; 
        ylabel("Imag($\sigma$)",'Interpreter','latex')
        %plot(ev_reps,evs(ev_reps),'o','MarkerFaceColor', 'b');
        %legend('approx.','single','Location','north')
        sgtitle(join([curve_param,', $\mu=$ ',num2str(mu),', \# of points = ',num2str(2*N+1)],""),'Interpreter','latex')
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        xlabel('N',"Interpreter","latex")
    end
    
    filename = ['figs/',curve_name,'_mu',num2str(mu),'_N',num2str(N),'.png'];
    filename = join(filename,'');
    %saveas(gcf,filename);
 
    % plot_ind = 1; %ev_reps(1);
    % figure()
    % plot3(real(evecs(:,plot_ind)),imag(evecs(:,plot_ind)),t(1:end-1));
    % xlabel('Re(u(t))')
    % ylabel('Im(u(t))')
    % zlabel('t')
    % 
    % figure()
    % plot(t(1:end-1),real(evecs(:,plot_ind)),'--o');
    % ylabel('Re(u(t))')
    % xlabel('t')
end    
% ----------------------------------------
%return 

%% Solution on M and in Omega ----------------------------------
fix = 11;
bessord = evsac_mat(evs_ac(fix)==evsac_mat(:,1),2); % for circle   
eden_fix = eden(:,fix);

%tic
    Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
    Rin = sqrt(Xindiff.^2+Yindiff.^2);
    HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
    Recon_integrand = HankRin.*transpose(eden_fix).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
    recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
    recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
    recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
    recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
    recint_sin = recint_sin1*sin_vecs(:,2:N)';
    recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
    %solution_fixed_point = 1i*pi/(4*N)*sum(Recon_integrand,2);
    solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
    %solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 
%toc


% %tic
% solution_fixed_point = zeros(max(size(points_inside)),1);
% for i = 1:length(solution_fixed_point)
%     diff1 = points_inside(:,i) - x(:,1:end-1);
%     diffR = mu*sqrt(diff1(1,:).^2+diff1(2,:).^2);
%     bessdiffR = besselh(0,diffR); bessdiffR = transpose(bessdiffR);
%     solution_fixed_point(i) = 1i*pi/(4*N)*sum(bessdiffR.*eden_fix.*dp2norm(fix,(1:end-1))');
% end
% %toc
% figure()
% subplot(1,2,1); plot3(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),'o')
% xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
% subplot(1,2,2); plot3(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),'or')
% xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Imag$(u)$','Interpreter','latex');

xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));
%[thetalin,rholin] = cart2pol(xlin,ylin);

[Xin2,Yin2] = meshgrid(xlin,ylin);
%[Rr, theta] = meshgrid(thetalin,rholin);

Uin = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
RealUin = real(Uin);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
%RealUin = RealUin/max(abs(solution_fixed_point));
ImagUin = imag(Uin);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');
%ImagUin = ImagUin/max(abs(solution_fixed_point));
bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;

% for circle only 
if curve_name == "circle"
    Rin2 = abs(mu)*sqrt(Xin2.^2+Yin2.^2); 
    BessRin2 = besselj(bessord,Rin2);
    cos_in = cos(bessord*atan2(Yin2,Xin2));
    sin_in = sin(bessord*atan2(Yin2,Xin2));
    Z_real_ac = BessRin2.*cos_in;
    Z_imag_ac = BessRin2.*sin_in;
    Z1 = Z_real_ac+1i*Z_imag_ac; Z2 = Z_real_ac-1i*Z_imag_ac; 

    app_sol_vec = reshape(Uin,[],1); 
    ac_sol1 = reshape(Z1,[],1); ac_sol1(isnan(app_sol_vec)) = [];
    ac_sol2 = reshape(Z2,[],1); ac_sol2(isnan(app_sol_vec)) = [];
    app_sol_vec(isnan(app_sol_vec)) = [];
    % ac_sol1 = reshape(Z_real_ac,[],1); ac_sol1(isnan(app_sol_vec)) = [];
    % ac_sol2 = reshape(Z_imag_ac,[],1); ac_sol2(isnan(app_sol_vec)) = [];
    %appacmat = [ac_sol1 ac_sol2 app_sol_vec];
    appacmat = [real(ac_sol1) imag(ac_sol1) real(ac_sol2) imag(ac_sol2) app_sol_vec];
    rrefappacmat = rref(appacmat);
    % rrefappacmat(1:5,:)
    ang1a = subspace(real(app_sol_vec), [real(ac_sol1) imag(ac_sol1)]);
    ang1ad = ang1a*180/pi;
    ang1b = subspace(imag(app_sol_vec), [real(ac_sol1) imag(ac_sol1)]);
    ang1bd = ang1b*180/pi;

    figure() 
    subplot(2,2,1);surf(Xin2,Yin2,RealUin);
    shading interp; view(0,90); colormap(magma); colorbar
    title('$Real(u^N_k)$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
    das = daspect; das(1) = das(2);
    set(gca,'DataAspectRatio',das);    

    subplot(2,2,2);surf(Xin2,Yin2,Z_real_ac);
    shading interp; view(0,90); colormap(magma); colorbar
    title('$Real(u_k)$','Interpreter','latex')
        das = daspect; das(1) = das(2);
    set(gca,'DataAspectRatio',das);    
    
    subplot(2,2,3);surf(Xin2,Yin2,ImagUin);
    shading interp; view(0,90); colormap(magma); colorbar
    title('$Imag(u^N_k)$','Interpreter','latex')
    xlabel('$x_1$','Interpreter','latex'); 
    ylabel('$x_2$','Interpreter','latex'); zlabel('Imag$(u)$','Interpreter','latex');
    das = daspect; das(1) = das(2);
    set(gca,'DataAspectRatio',das);    

    subplot(2,2,4);surf(Xin2,Yin2,Z_imag_ac);
    shading interp; view(0,90); colormap(magma); colorbar
    xlabel('$x_1$','Interpreter','latex'); 
    title('$Imag(u_k)$','Interpreter','latex')
    das = daspect; das(1) = das(2);
    set(gca,'DataAspectRatio',das);    
    
    sgtitle({join([curve_name,', $\mu=$ ',num2str(mu),', N = ',num2str(N),', k = ',num2str(fix),', $\sigma_k$ = ',num2str(evs(fix))]) ...
        join(['Angle between $Real (u^N_k)$ and $Real,Imag (u_k)$ = ',num2str(ang1a),' rad = ',num2str(ang1ad),' deg']) ...
        join(['Angle between $Imag (u^N_k)$ and $Real,Imag (u_k)$ = ',num2str(ang1b),' rad = ',num2str(ang1bd),' deg'])
        },'Interpreter','latex','fontsize',18)

    disp(['angle between real part of approx sol and real and imag parts of actual sol = ',num2str(ang1a),' rad = ',num2str(ang1ad),' deg']); 
    set(findobj(gcf,'type','axes'),'FontSize',18)
    filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/efn_disk/',curve_name,num2str(disk_r),'_mu',num2str(mu),'_N_',num2str(N),'_k',num2str(fix),'_m',num2str(m),'_inside.png'],"")];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 25 25]); 
    %exportgraphics(gcf,filename,'Resolution',300)

    % app_sol_vec = reshape(RealUin,[],1); 
    % ac_sol1 = reshape(Z_real_ac,[],1); ac_sol1(isnan(app_sol_vec)) = [];
    % app_sol_vec(isnan(app_sol_vec)) = [];
    % appacmat = [ac_sol1 app_sol_vec];
    % rrefappacmat = rref(appacmat);
    % rrefappacmat(1:5,:)

    ang2a = subspace(real(evecs(:,fix)),[cos(bessord*t(1:end-1))' sin(bessord*t(1:end-1))']);
    ang2ad = ang2a*180/pi;
    ang2b = subspace(imag(evecs(:,fix)),[cos(bessord*t(1:end-1))' sin(bessord*t(1:end-1))']);
    ang2bd = ang2b*180/pi;

    figure() 
    subplot(2,2,1)
    plot(t(1:end-1),real(eden(:,fix)),'LineWidth',3);
   % ylim([-1 1])
    %title('Density','Interpreter','latex')
    xlabel('$t$','Interpreter','latex'); ylabel('$Real(\psi^N_k)$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
    
    subplot(2,2,2); 
    plot(t(1:end-1),real(evecs(:,fix)),'r','LineWidth',3);
   % ylim([-1 1])
    %title('Eigen function','Interpreter','latex')
    xlabel('$t$','Interpreter','latex'); ylabel('$Real(u^N_k)$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
    
    subplot(2,2,3)
    plot(t(1:end-1),imag(eden(:,fix)),'LineWidth',3);
    %ylim([-1 1])
    %title('Density','Interpreter','latex')
    xlabel('$t$','Interpreter','latex'); ylabel('$Imag(\psi^N_k)$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');
    
    subplot(2,2,4); 
    plot(t(1:end-1),imag(evecs(:,fix)),'r','LineWidth',3);
    %ylim([-1 1])
    %title('Eigen function','Interpreter','latex')
    xlabel('$t$','Interpreter','latex'); ylabel('$Imag(u^N_k)$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');

    sgtitle({join(['On the boundary of ',curve_name,' r = ',len/(2*pi),', $\mu=$ ',num2str(mu),', N = ',num2str(N),', k = ',num2str(fix),', $\sigma_k$ = ',num2str(evs(fix))]) ...
        join(['Angle between $Real (u^N_k)$ and $Real,Imag (u_k)$ = ',num2str(ang2a),' rad = ',num2str(ang2ad),' deg']) ...
        join(['Angle between $Imag (u^N_k)$ and $Real,Imag (u_k)$ = ',num2str(ang2b),' rad = ',num2str(ang2bd),' deg'])
        },'Interpreter','latex','fontsize',18)
      set(findobj(gcf,'type','axes'),'FontSize',18)
    filename = [join(['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/efn_disk/',curve_name,num2str(disk_r),'_mu',num2str(mu),'_N_',num2str(N),'_k',num2str(fix),'_m',num2str(m),'_onM.png'],"")];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 25 25]); 
    %exportgraphics(gcf,filename,'Resolution',300)
else
    figure() 
    subplot(3,2,[1,3]);surf(Xin2,Yin2,RealUin);
    shading interp; view(0,90); colormap(magma); colorbar
    title('Real(U)','Interpreter','latex')
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Real$(u)$','Interpreter','latex');   
    
    subplot(3,2,[2,4]);surf(Xin2,Yin2,ImagUin);
    shading interp; view(0,90); colormap(magma); colorbar
    title('Imag(U)','Interpreter','latex')
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); zlabel('Imag$(u)$','Interpreter','latex'); 
    sgtitle(join([curve_param,', $\mu=$ ',num2str(mu),', N = ',num2str(N),', k = ',num2str(fix),', $\sigma_k$ = ',num2str(evs(fix))]),'Interpreter','latex')
    
end

% Densities ------------
% figure()
% %plot(t(1:end-1),real(evecs(:,plot_ind)),'r',t(1:end-1),imag(evecs(:,plot_ind)),'b'); hold on;
% plot(t(1:end-1),real(rescale(real(evecs(:,plot_ind)),-1,1)),'r',t(1:end-1),imag(rescale(real(evecs(:,1)),-1,1)),'b'); hold on;
% plot(t(1:end-1),cos(bessord*t(1:end-1)),'c',t(1:end-1),sin(bessord*t(1:end-1)),'k');
% xlabel('t')
% ylabel('u')
% title('lambda=',num2str(evs(plot_ind)))
    
figure()
subplot(2,1,1);
plot3(t(1:end-1),real(eden(:,fix)),imag(eden(:,fix))); hold on;
xlabel('t')
ylabel('R(u)')
zlabel('Im(u)')
title('density')

subplot(2,1,2);
plot3(t(1:end-1),real(evecs(:,fix)),imag(evecs(:,fix))); hold on;

if curve_name == "circle"
plot3(t(1:end-1),cos(bessord*t(1:end-1)),sin(bessord*t(1:end-1)));
end

ang2 = subspace((evecs(:,fix)),[cos(bessord*t(1:end-1))' sin(bessord*t(1:end-1))']);
ang2d = ang2*180/pi;
disp(['angle between approx efn on M and real and imag parts of actual efn on M = ',num2str(ang2),' rad = ',num2str(ang2d),' deg']); 

xlabel('t')
ylabel('R(u)')
zlabel('Im(u)')
title('efn on M')
sgtitle(['sigma=',num2str(evs(fix))])

% figure()
% subplot(1,2,1); plot3(x(1,1:end-1),x(2,1:end-1),real(evecs(:,fix)))
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% zlabel('$Re(u)$','Interpreter','latex')
% title('imag part of efn on M')
% 
% subplot(1,2,2); plot3(x(1,1:end-1),x(2,1:end-1),imag(evecs(:,fix)))
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% zlabel('$Im(u)$','Interpreter','latex')
% title('imag part of efn on M')

figure()
%plot3(x(1,1:end-1),x(2,1:end-1),real(evecs(:,fix))); hold on;
%plot3(x(1,1:end-1),x(2,1:end-1),imag(evecs(:,fix)))
plot3(x(1,:),x(2,:),[real(evecs(:,fix));real(evecs(1,fix))]); hold on;
plot3(x(1,:),x(2,:),[imag(evecs(:,fix));imag(evecs(1,fix))]); hold on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
zlabel('$u$','Interpreter','latex')
title(['$\sigma$',num2str(fix),' = ',num2str(evs(fix))],'Interpreter','latex')
legend(["Real", "Imag"],'Location','north')

evs(1) < (2-mu^2)*pi/(len)

%% Weyl --------------------------------
% asymp_evs = uniquetol(evs);
% asymp_evs = asymp_evs(1:N/2);
asymp_evs = (evs(1:N)); asymp_evs1 = asymp_evs;
nN = length(asymp_evs);

start = 1; 
% if mu > 2
% stop = round(mu^2);
% else
stop = N;round(mu*pi/dom_area); round(mu^2*pi/(2*pi*dom_area));
% end    
%stop = 500;%length(asymp_evs);
y = asymp_evs(start:stop)*len;
x1 = start:stop;
%w = [10;zeros(length(x1)-2,1)+1;10];
%w = [10;zeros(length(x1)-1,1)+1];

% myfit = fittype('a + b*log(x)+c*exp(d*x)',...
% 'dependent',{'y'},'independent',{'x'},...
% 'coefficients',{'a','b','c','d'});
% 
% logsf = fit(x1',y,myfit);
% %logsf = fit(x1',y,myfit,'Weight',w);
% 
% figure();
% plot(start:stop,asymp_evs(start:stop)*len); hold on;
% plot(start:stop,logsf(start:stop)); hold on;

myfit = fittype('a+b*(x)+c/sqrt(x)+d*sqrt(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b','c','d'});

% myfit1 = fittype('a + b*sqrt(x)+c*log(x)',...
% 'dependent',{'y'},'independent',{'x'},...
% 'coefficients',{'a','b','c'});
% 
% myfit2 = fittype('a + b*x',...
% 'dependent',{'y'},'independent',{'x'},...
% 'coefficients',{'a','b'});
% 
% myfit3 = fittype('a + b*log(x)',...
% 'dependent',{'y'},'independent',{'x'},...
% 'coefficients',{'a','b'});
% 
% myfit4 = fittype('a + b*sqrt(x)',...
% 'dependent',{'y'},'independent',{'x'},...
% 'coefficients',{'a','b'});
% 
% myfit5 = fittype('a*x + b*sqrt(x)',...
% 'dependent',{'y'},'independent',{'x'},...
% 'coefficients',{'a','b'});

expsf = fit(x1',y,'exp1');
logsf = fit(x1',y,myfit);
% logsf1 = fit(x1',y,myfit1);
% logsf2 = fit(x1',y,myfit2);
% logsf3 = fit(x1',y,myfit3);
% logsf4 = fit(x1',y,myfit4);
% logsf5 = fit(x1',y,myfit5);
figure()
subplot(3,1,1)
plot(start:stop,asymp_evs(start:stop)*len,'--*','LineWidth',2); hold on;
ylabel("$|M|\sigma_k$",'Interpreter','latex')
xlabel("$k$",'Interpreter','latex')
%plot(1:nN,asymp_evs(1:nN)*len,'--*','LineWidth',2); hold on;

% plot(start:stop,expsf(start:stop),'LineWidth',1.6); hold on;
plot(start:stop,logsf(start:stop),'LineWidth',3); hold on;
legend(["$|M|\sigma_k$",join([logsf.c,"$/\sqrt k+$",logsf.b,"$k$"],"")],'Interpreter','latex')
% plot(start:stop,logsf1(start:stop),'LineWidth',1.6); hold on;
% plot(start:stop,logsf2(start:stop),'LineWidth',1.6); hold on;
% plot(start:stop,logsf3(start:stop),'LineWidth',1.6); hold on;
% plot(start:stop,logsf4(start:stop),'LineWidth',1.6); hold on;
% plot(start:stop,logsf5(start:stop),'LineWidth',1.6); hold on;
% legtext = ["$truth$","$e^x$","$\log x+x^2$","$\log x+\sqrt{x}$","line","$\log x$","$\sqrt{x}$","$x+\sqrt{x}$"];
% legend(legtext,'interpreter','latex')

l1 = polyfit(stop:nN,(asymp_evs(stop:nN))*len,1); l1
l1 = real(l1);
% figure(); 
subplot(3,1,2)
plot(stop:nN,(asymp_evs(stop:nN))*len,'--*','LineWidth',2); hold on;
plot(stop:nN,l1(1)*(stop:nN)+l1(2),'LineWidth',3); 
xlabel("$k$",'Interpreter','latex')
ylabel("$|M|\sigma_k$",'Interpreter','latex')
legend(["$|M|\sigma_k$",join([l1(2),"$+$",l1(1),"$k$"],"")],'Interpreter','latex')
%legend(["$|M|\sigma_k$",join([logsf.c,"$/\sqrt k+$",logsf.b,"$k$"],""),join([l1(2),"$+$",l1(1),"$k$"],"")],'Interpreter','latex')
%title(join([curve_name,", $|M| = $ ",len,", $\mu = $ ",mu],""),'Interpreter','latex')

sgtitle(join([curve_name,", $|M| = $ ",len,", $\mu = $ ",mu],""),'Interpreter','latex','FontSize',20)

set(findobj(gcf,'type','axes'),'FontSize',20)
% logsf.b/mu^2
% mu^2/logsf.b

mu/(pi/2)

asymp_evs = sort(abs(asymp_evs));
%asymp_evs = uniquetol(asymp_evs);
lambda = 1:N;
N_Lambda = zeros(length(N),1);
for nl = 1:N%length(asymp_evs)
    %N_Lambda(nl) = sum(asymp_evs*len < nl);
    N_Lambda(nl) = sum(asymp_evs(1:nl)*len <= nl);
end    

myfit_w = fittype('b*(x)+c*sqrt(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'b','c'});

myfit_w1 = fittype('b*sqrt(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'b'});

% dees = diff(N_Lambda);
% dees = find(dees < 1); dees = dees(dees > 100);
% if mu ~= 1
%     Stop = dees(1);
% else
%    Stop = length(N_Lambda);
% end

Stop = stop;%mu^2;%length(N_Lambda);
% end

Start = 1; 

weyl_c = fit((Start:N)',N_Lambda(Start:N)',myfit_w);
weyl_c1 = fit((Start:Stop)',N_Lambda(Start:Stop)',myfit_w);
% weyl_c2 = fit((Stop:N)',N_Lambda(Stop:N)',myfit_w);
% weyl_2 = weyl_c1.b;

weyl_have = [weyl_c.b,weyl_c.c];
%weyl_want = [dom_area/(pi),len/(4*pi)];
weyl_want = [len/(pi),len/(4*pi)];

% weyl_have./weyl_want
% weyl_want./weyl_have
% weyl_want./weyl_have*mu
% weyl_want./weyl_have/mu
% weyl_want./weyl_have/mu^2


%figure(); 
% subplot(2,2,1)
% plot(Start:Stop,N_Lambda(Start:Stop),'--o','LineWidth',2); hold on;
% plot(Start:Stop,weyl_c1(Start:Stop),'--*','LineWidth',1); hold on;
% 
% subplot(2,2,2)
% plot(Stop:N,N_Lambda(Stop:N),'--o','LineWidth',2); hold on;
% plot(Stop:N,weyl_c2(Stop:N),'--*','LineWidth',1); hold on;
% 
% subplot(2,2,3:4)
subplot(3,1,3)
plot(Start:N,N_Lambda(Start:N),'--o','LineWidth',2); hold on;
plot(Start:N,weyl_c(Start:N),'--*','LineWidth',1); hold on;
ylabel("$N(|M||\sigma_s|)$",'Interpreter','latex')
xlabel("$s$",'Interpreter','latex')
legend(["$N(|M||\sigma_s|)$",join([weyl_c.c,"$\sqrt s+$",weyl_c.b,"$s$"],"")],'Interpreter','latex')%plot(Start:Stop,weyl_want(1)*(Start:Stop)+weyl_want(2)*sqrt(Start:Stop),'LineWidth',2); hold on;
title(join(["$1/$",weyl_c.b," $\approx$ ",1/weyl_c.b],""),'Interpreter','latex')

set(findobj(gcf,'type','axes'),'FontSize',20)

filename = join(['august talk/weyl/',curve_name,'_mu_',num2str(mu),'_N_',N,'.png'],"");
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
ax = gcf;
%exportgraphics(ax,filename,'Resolution',300) 
%savefig(ax,strrep(filename,'.png','.fig'))

disp("Weyl have = "+weyl_have)
disp("2area/(mu*pi) = " +1/(mu/dom_area*pi/2))
weyl_have(2)/(dom_area*mu/pi^2)
% weyl_c2.c*mu^2
% weyl_want/weyl_c2.c

% % plot(start:stop,logsf1(start:stop)); hold on;
% plot(start:stop,expsf(start:stop));

% %disp(mu^2*2*pi); 
% disp(logsf); 
% 
% log_coeffs = [ones(size(x1')),log(x1')]\y;
% disp(log_coeffs(2)/(stop))
%toc

%% efun normalization on M ----------------------
fix = 1;
norm_M = sum(abs(evecs(:,fix)).^2)*pi/N;
norm_evec = evecs(:,fix)/sqrt(abs(norm_M));
% figure();
% subplot(2,1,1); plot(real(norm_evec))
% subplot(2,1,2); plot(imag(norm_evec))
norm_M1 = sum((norm_evec).^2)*pi/N;
bessord = evsac_mat(evs_ac(fix)==evsac_mat(:,1),2); % for circle   
eden_fix = eden(:,fix)/sqrt(abs(norm_M));

norm_evec1 = B*eden_fix/2;
norm_normder_evecs = (A+I)*eden_fix/2;

% point_ind = 10;
% point = x(:,point_ind); nx(:,point_ind)
% h = 1e-3; in_point_x = point - h*[1;0]; 
% in_point_y = point - h*[0;1]; 
% 
% for i = 1:2
% % these have the same names as in computing solution inside everywhere, sorrs
% if i == 1
% [Xtauin,Xin] = meshgrid(x(1,1:end-1),in_point_x(1,:)); % mesh x coord of in points and bdry points
% [Ytauin,Yin] = meshgrid(x(2,1:end-1),in_point_x(2,:)); % mesh y coord of in points and bdry points
% else
% [Xtauin,Xin] = meshgrid(x(1,1:end-1),in_point_y(1,:)); % mesh x coord of in points and bdry points
% [Ytauin,Yin] = meshgrid(x(2,1:end-1),in_point_y(2,:)); % mesh y coord of in points and bdry points
% end    
% Xindiff = Xin-Xtauin; Yindiff = Yin-Ytauin;
% Rin = sqrt(Xindiff.^2+Yindiff.^2);
% HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
% Recon_integrand = HankRin.*eden_fix'.*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
% recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
% recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
% recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
% recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
% recint_sin = recint_sin1*sin_vecs(:,2:N)';
% recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
% if i==1
% sol_inpointx = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
% else
% sol_inpointy = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
% end
% end
% sol_point = norm_evec(point_ind);
% jump = eden_fix(point_ind)/2;
% [sol_point ,sol_inpointx,sol_inpointy]
% app_grad = [(sol_point-sol_inpointx)/h,(sol_point-sol_inpointy)/h]; 
% 
% [sum(app_grad.*nx(:,point_ind)'), norm_normder_evecs(point_ind)]
% 
% % if we normalize with taking abs of norm_M, then
% % the imag part goes away !
% 
% figure(); 
% plot(x(1,:),x(2,:),'--o'); axis equal
% hold on;
% 
% exT = max(max(abs(x)));
% lenT = 2*exT;
% grid1 = -lenT:(2*lenT)/99:lenT;
% [GridY,GridX] = meshgrid(grid1,grid1);
% plot(GridX,GridY,'*'); hold on;
% 
% idx1 = inpolygon(GridX(:),GridY(:),x(1,:),x(2,:)) ;
% points_inside1 = transpose([GridX(idx1),GridY(idx1)]);
% plot(points_inside1(1,:),points_inside1(2,:),'o')
% 
% N_omega = size(points_inside1,2);
% n_lambda = ceil(sqrt(N_omega/4));
% if rem(n_lambda,2) ~=1
%     n_lambda = n_lambda-1;
% end    
% finds = (n_lambda-1)/2; finds = -finds:finds;
% [Finds2,Finds1] = meshgrid(finds);  % Phi_j = exp(i x\dot I 2pi/2T)
% % x\dot I = dot product of point with fourier index
% 
% N_lambda = n_lambda^2; 
% N_omega/N_lambda % what happens if < 4 ?
% 
% Fext_A = zeros(N_omega,N_lambda);
% 
% const = 1i*2*pi/(lenT);
% for k = 1:(N_omega) 
%     %for j = 1:length(N_lambda)
%         Fext_A(k,:) = exp(const*reshape((points_inside1(1,k)*Finds1+points_inside1(2,k)*Finds2)',1,N_lambda)); 
%     %end
% end    
% Fext_A = Fext_A/(2*lenT*length(grid1)); %
% 
% [Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside1(1,:)); % mesh x coord of in points and bdry points
% [Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside1(2,:)); % mesh y coord of in points and bdry points
% 
% %tic
%     Xindiff = Xin-Xtauin; Yindiff = Yin-Ytauin;
%     Rin = sqrt(Xindiff.^2+Yindiff.^2);
%     HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
%     Recon_integrand = HankRin.*eden_fix'.*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
%     recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
%     recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
%     recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
%     recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
%     recint_sin = recint_sin1*sin_vecs(:,2:N)';
%     recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
%     %solution_fixed_point = 1i*pi/(4*N)*sum(Recon_integrand,2);
%     solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
%     %solution_fixed_point_mat = 1i*pi/(4*N)*arrayfun(@(x) trapz([Recon_integrand(x,:),Recon_integrand(x,1)]), 1:size(recint,1)); 
% %toc
% 
% P = Fext_A*(Fext_A)'-eye(N_omega);
% a1 = (P*Fext_A)\(P*solution_fixed_point);
% a2 = (Fext_A)'*(solution_fixed_point-Fext_A*a1);
% Fourier_coeffs = a1+a2;
% 
% gradF = @(point_nmbr) const/length(grid1)*[sum(exp(const*reshape((x(1,point_nmbr).*Finds1+x(2,point_nmbr).*Finds2)',1,N_lambda)).*Fourier_coeffs'.*reshape(Finds1',1,N_lambda)),sum(exp(const*reshape((x(1,point_nmbr).*Finds1+x(2,point_nmbr).*Finds2)',1,N_lambda)).*Fourier_coeffs'.*reshape(Finds2',1,N_lambda))];
% gradF(2)
% 
% tang_Grad = @(point_nmbr) gradF(point_nmbr)-transpose(sum(gradF(point_nmbr).*nx(:,point_nmbr)')*nx(:,point_nmbr));
% 
% nabla_T = zeros(size(x,2),2); nabla = zeros(size(x,2),2);
% for i = 1:size(x,2)
% nabla_T(i,:) = tang_Grad(i); 
% nabla(i,:) = gradF(i);
% end
% vecnorm(transpose(nabla_T));
% subspace([real(tan_grad_disk) imag(tan_grad_disk)],imag(nabla_T))

tx = [-nx(2,:);nx(1,:)]; %tx.*nx
tan_grad_disk = transpose(1i/disk_r*bessord*besselj(bessord,mu*disk_r)*exp(1i*bessord*t)).*tx';
%tan_grad_disk = (1i/disk_r*bessord*besselj(bessord,mu*disk_r)*exp(1i*bessord*t))'.*tx';
% [tan_grad_disk./vecnorm(transpose(tan_grad_disk))', nabla_T./vecnorm(transpose(nabla_T))']
% [vecnorm(transpose(tan_grad_disk));vecnorm(transpose(nabla_T))]'

tan_grad_disk = tan_grad_disk(1:end-1,:); 

%% rbf extension ------------
% figure(); 
% plot(x(1,:),x(2,:),'--o'); axis equal
% hold on;

% d1 = 0.1; d2 = 1;  
% points_inside_close = x(:,1:end-1) - d1*nx(:,1:end-1);
% points_outside_close = x(:,1:end-1) + d2*nx(:,1:end-1);
% points_inside_close = x(:,1:2:end-1) - d1*nx(:,1:2:end-1);
% points_outside_close = x(:,1:2:end-1) + d2*nx(:,1:2:end-1);
% plot(points_inside_close(1,:),points_inside_close(2,:),'--o'); axis equal
% hold on;
% plot(points_outside_close(1,:),points_outside_close(2,:),'--o'); axis equal

% rbf_points = [points_inside_close,x(:,1:end-1),points_outside_close];
% rbf_f = [norm_evec-d1*(nx(1,1:end-1)'+1i*nx(2,1:end-1)');norm_evec;norm_evec+d2*(nx(1,1:end-1)'+1i*nx(2,1:end-1)')];
% rbf_points = [x(:,1:end-1),points_outside_close];
%rbf_f = [norm_evec;norm_evec+d2*(nx(1,1:end-1)'+1i*nx(2,1:end-1)')];
rbf_points = x(:,1:end-1);
rbf_f = norm_evec;
%rbf_f = transpose(exp(1i*bessord*t(1:end-1)))*besselj(bessord,mu*disk_r);
%rbf_f = [real(norm_evec)-d1*(nx(2,1:end-1)');imag(norm_evec)+d2*(nx(2,1:end-1)')];
%rbf_f = [real(norm_evec)-d1*(nx(1,1:end-1)');norm_evec+d2*(nx(1,1:end-1)')];
% rbf_points = [points_outside_close];
% rbf_f = norm_evec+d2*(nx(1,1:end-1)'+1i*nx(2,1:end-1)');
%rbf_f = [norm_evec(1:2:end)-d1*(nx(1,1:2:end-1)'+1i*nx(2,1:2:end-1)');norm_evec(1:2:end)+d2*(nx(1,1:2:end-1)'+1i*nx(2,1:2:end-1)')];
%rbf_f = [norm_evec-d1*(1+1i);norm_evec+d2*(1+1i)];

[Xtau_rbf,Xt_rbf] = meshgrid(rbf_points(1,:));
[Ytau_rbf,Yt_rbf] = meshgrid(rbf_points(2,:));
R_rbf = sqrt((Xt_rbf-Xtau_rbf).^2+(Yt_rbf-Ytau_rbf).^2);

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

p_deg = 3; p_basis_size = p_deg*(p_deg+1)/2+p_deg+1;
P = zeros(size(rbf_points,2),p_basis_size);
for i = 1:size(rbf_points,2)
    if p_deg == 1
        P(i,:) = [1;rbf_points(1,i);rbf_points(2,i)];
    elseif p_deg == 2
        P(i,:) = [1;rbf_points(1,i);rbf_points(2,i);rbf_points(1,i)*rbf_points(2,i);rbf_points(1,i)^2;rbf_points(2,i)^2];
    elseif p_deg == 3
        P(i,:) = [1;rbf_points(1,i);rbf_points(2,i);...
            rbf_points(1,i)*rbf_points(2,i);rbf_points(1,i)^2;rbf_points(2,i)^2;...
            rbf_points(1,i)^2*rbf_points(2,i);rbf_points(1,i)*rbf_points(2,i)^2;...
            rbf_points(1,i)^3;rbf_points(2,i)^3
            ];
    end    
end    

% rbf_mat = zeros(size(rbf_points,2)+p_basis_size);
% rbf_mat(1:size(rbf_points,2),1:size(rbf_points,2)) = A_rbf;
% rbf_mat(1:size(rbf_points,2),(size(rbf_points,2)+1):end) = P;
% rbf_mat((size(rbf_points,2)+1):end,1:size(rbf_points,2)) = P';
% rbf_mat = rbf_mat/sqrt(size(rbf_points,2));
rbf_mat = A_rbf;

%rbf_rhs = [rbf_f;zeros(p_basis_size,1)];
rbf_rhs = rbf_f;
%rbf_coeffs = rbf_mat\rbf_lhs;
rbf_coeffs = lsqminnorm(rbf_mat,rbf_rhs);
%poly_coefs = rbf_coeffs((end-p_basis_size+1):end);

rbf_grad = zeros(size(x,2)-1,2);
rbf_fun = zeros(size(x,2)-1,1);
rbf_tgrad = zeros(size(x,2)-1,2);

for i = 1:(size(x,2)-1)
    for j = 1:size(rbf_points,2)
        if bf == 1
            rbf_fun(i) = rbf_fun(i)+rbf_coeffs(j)*rbf(x(1,i)-rbf_points(1,j),x(2,i)-rbf_points(2,j));
            rbf_grad(i,:) = rbf_grad(i,:) + rbf_coeffs(j)*[drbfx(x(1,i)-rbf_points(1,j),x(2,i)-rbf_points(2,j)),drbfy(x(1,i)-rbf_points(1,j),x(2,i)-rbf_points(2,j))]; 
        else
            distr = norm(x(:,i)-rbf_points(:,j)); 
            rbf_fun(i) = rbf_fun(i)+rbf_coeffs(j)*rbf(distr);
            rbf_grad(i,:) = rbf_grad(i,:) + rbf_coeffs(j)*drbf(distr)*(x(:,i)'-rbf_points(:,j)'); 
        end


        % if p_deg == 1
        %     rbf_grad(i,:) = rbf_grad(i,:) + [poly_coefs(2),poly_coefs(3)];
        % elseif p_deg == 2
        %     rbf_grad(i,:) = rbf_grad(i,:) + [poly_coefs(2)+poly_coefs(4)*x(2,i)+2*poly_coefs(5)*x(1,i),poly_coefs(3)+poly_coefs(4)*x(1,i)+2*poly_coefs(6)*x(2,i)];
        % elseif p_deg == 3
        %     rbf_grad(i,:) = rbf_grad(i,:) + [poly_coefs(2)+poly_coefs(4)*x(2,i)+2*poly_coefs(5)*x(1,i)+2*poly_coefs(7)*x(1,i)*x(2,i)+poly_coefs(8)*x(2,i)^2+poly_coefs(9)*x(1,i)^3,...
        %         poly_coefs(3)+poly_coefs(4)*x(1,i)+2*poly_coefs(6)*x(2,i)+poly_coefs(7)*x(1,i)^2+2*poly_coefs(8)*x(1,i)*x(2,i)+3*poly_coefs(10)*x(2,i)^2];
        % end    
    end    
    rbf_tgrad(i,:) = rbf_grad(i,:)-sum(rbf_grad(i,:).*nx(:,i)')*nx(:,i)';
end

% blahblah2 = [tan_grad_disk./vecnorm(transpose(tan_grad_disk))', rbf_tgrad./vecnorm(transpose(rbf_tgrad))']
% blahblah = [tan_grad_disk, rbf_tgrad];
% 
% sum(abs(rbf_fun-norm_evec)<1e-12)
% figure()
% subplot(1,2,1); plot(real(tan_grad_disk))
% subplot(1,2,2); plot(real(rbf_tgrad))
% figure(); subplot(1,2,1); plot(imag(tan_grad_disk))
% subplot(1,2,2); plot(imag(rbf_tgrad))
% 
% figure()
% subplot(2,2,1)
% plot(imag(norm_evec))
% hold on; plot(sum(imag(tan_grad_disk).*tx(:,1:end-1)',2))
% subplot(2,2,2)
% plot(imag(rbf_fun))
% hold on; plot(sum(imag(rbf_tgrad).*tx(:,1:end-1)',2))
% 
% %figure()
% subplot(2,2,3)
% plot(real(norm_evec))
% hold on; plot(sum(real(tan_grad_disk).*tx(:,1:end-1)',2))
% subplot(2,2,4)
% plot(real(rbf_fun))
% hold on; plot(sum(real(rbf_t  grad).*tx(:,1:end-1)',2))

subspace(real(rbf_tgrad),[real(tan_grad_disk) imag(tan_grad_disk)])
subspace(imag(rbf_tgrad),[real(tan_grad_disk) imag(tan_grad_disk)])

v = [cos(t);sin(t)]; v = v(:,1:end-1);  
vn = sum(v.*nx(:,1:end-1)); dz = dp2norm(fix,(1:end-1));
dndz = vn.*dz;
tang_normsq = vecnorm(transpose(rbf_tgrad)).^2;
integrand1 = tang_normsq.*dndz;
int1 = sum(integrand1)*pi/N;
evk = evs(fix);
integrand2 = (mu^2+evk^2+kappa*evk)*abs(norm_evec).^2; 
int2 = sum(integrand2)*pi/N;

Int1 = len*(int1-int2); 
Int2 = evk*kappa*sum(dndz)*pi/N; 
Int = Int1+Int2;

%% mu --> Laplace DEV ---
[U,S,V] = svd(B); %toc
S = diag(S); [s_min,min_ind] = min(S); 
%all_S(:,ms) = S;
fix = min_ind;
phi_min = V(:,fix);
Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
Rin = sqrt(Xindiff.^2+Yindiff.^2);
HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
Recon_integrand = HankRin.*transpose(phi_min).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); 
%mat_trapz = trapz(t,[recint,recint(:,1)],2)*1i/4

xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));
%[thetalin,rholin] = cart2pol(xlin,ylin);

[Xin2,Yin2] = meshgrid(xlin,ylin);
%[Rr, theta] = meshgrid(thetalin,rholin);

Win = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;
Win(~bdrypts) = NaN;

figure(); plot3(t(1:end-1),real(B*phi_min/2),imag(B*phi_min/2),'*')

k = 1; 
evec_k = evecs(:,k);
eden_fix = eden(:,k);

Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
Rin = sqrt(Xindiff.^2+Yindiff.^2);
HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
Recon_integrand = HankRin.*transpose(eden_fix).*dp2norm(k,(1:end-1)); % make this a tp and then sum and scale
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); 
%mat_trapz = trapz(t,[recint,recint(:,1)],2)*1i/4

xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));
%[thetalin,rholin] = cart2pol(xlin,ylin);

[Xin2,Yin2] = meshgrid(xlin,ylin);
%[Rr, theta] = meshgrid(thetalin,rholin);

Uin = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;
Uin(~bdrypts) = NaN;
RealUin = real(Uin);%griddata(points_inside(1,:),points_inside(2,:),real(solution_fixed_point),Xin2,Yin2,'cubic');
ImagUin = imag(Uin);%griddata(points_inside(1,:),points_inside(2,:),imag(solution_fixed_point),Xin2,Yin2,'cubic');

Uin_magsq = abs(Uin).^2;
Uin_magsq(isnan(Uin_magsq)) = 0;
I = trapz(ylin,trapz(xlin,Uin_magsq,2));
Uin = Uin/I;

UW = Uin.*Win;
UW(isnan(UW)) = 0;
I2 = trapz(ylin,trapz(xlin,UW,2));

toc