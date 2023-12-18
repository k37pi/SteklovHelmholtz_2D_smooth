clear; close all; 
%% inputs -----------------------------------------------
mu_min = 7.5; mu_max = 10.5; 
curve_param = "egg0.2";
eps = 0.2; N = 240; RI = 2;
step =  0.010033;0.0033445;0.0016722;0.011136;

filename = 'files/transmission/';%['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
%filename = join([filename,curve_param,"_ri_",RI,"min_phi_",mu_min,"_",mu_max,"_N",num2str(N)],'');
filename = join([filename,curve_param,"_ri_",RI,"min_phi_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".csv"],"");
phis = readmatrix(filename);

filename = 'files/transmission/';%['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
%filename = join([filename,curve_param,"_ri_",RI,"minsigmu_",mu_min,"_",mu_max,"_N",num2str(N)],'');
filename = join([filename,curve_param,"_ri_",RI,"minsigmu_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".csv"],"");
musig = readmatrix(filename);

filename = 'files/transmission/';%['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
%filename = join([filename,curve_param,"_ri_",RI,"min_inds_",mu_min,"_",mu_max,"_N",num2str(N)],'');
filename = join([filename,curve_param,"_ri_",RI,"min_inds_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".csv"],"");
mininds = readmatrix(filename);
 
mus = musig(musig(:,1)>=mu_min & musig(:,1)<=mu_max,1);
sigs = musig(musig(:,1)>=mu_min & musig(:,1)<=mu_max,2);
phis = phis(:,musig(:,1)>=mu_min & musig(:,1)<=mu_max);
min_inds = mininds(musig(:,1)>=mu_min & musig(:,1)<=mu_max);

mus = transpose(mus); sigs = transpose(sigs);
s_mins = sigs;

figure()
plot(mus,s_mins,'--o','LineWidth',2)
xlabel('$\mu$','Interpreter','latex')
ylabel('$s_{\min}$','Interpreter','latex')
title(join(["Minimum singular values of $C(\mu)$ for ",curve_param,", $n_{RI}$= ",RI,", N = ",num2str(N)],""),'Interpreter','latex')
filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,curve_param,"min_sig_",mu_min,"_",mu_max,"_N",num2str(N),"_step",num2str(step)],'');
filename = strrep(filename,'.',"dot");
filename = join([filename,".png"],"");
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
set(findobj(gcf,'type','axes'),'FontSize',18)


I = 1:length(mus); % all lam points
J = I(2:end-1); % interior points
S = sigs; 
J = J( S(J)<S(J-1) & S(J)<S(J+1) ); % local minima
J= J+ (S(J-1)>S(J+1)); % points where sign changes
K = 0*I; K(J) = 1;
S = S.*(-1).^cumsum(K); % introduce sign flips
hold off, %plot(mus,S,'*',"linewidth",2), hold on % plot signed angle function
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
x = uniquetol(x,1e-4);
x(isnan(x)) = [];

hold on; plot(points,rbf_fun,'color',"#D95319",'LineWidth',3.5)
plot(mus,s_mins,'-','color',"#0072BD",'LineWidth',2)
legend(["$s_{min}$", "RBF"],"AutoUpdate","off","Interpreter","latex")
xline(uniquetol(x),'k',"linewidth",1.5);
yline(0)
% % ---
% return
% Lamb = zeros(length(J),1);
% % Find eigenvalues via local 3rd-order interpolation works... why tf ?
% for j = 1:length(J)
% % I = J(j)-2:J(j)+1; % for this config .. ? ri = 4, 3<= mus <= 7
% % lambda = polyval(polyfit(S(I)/norm(S(I)),mus(I),3),0);
% I = J(j)-1:J(j)+1; % for this config .. ? ri = 4, 3<= mus <= 7
% lambda = polyval(polyfit(S(I)/norm(S(I)),mus(I),2),0);
% % I = J(j)-2:J(j)+1; % for this config .. ? ri = 4, 3<= mus <= 7
% % lambda = polyval(polyfit(S(I)/norm(S(I)),mus(I),3),0);
% disp(lambda) % display eigenvalue
% Lamb(j) = lambda;
% plot(lambda*[1 1],.8*[-1 1],'r') % plot eigenvalue
% end

xlim([min(mus) max(mus)*1.0])
xticks(uniquetol(x))
xtickangle(80)
%set(gcf, 'PaperUnits', 'centimeters');
%set(gcf, 'PaperPosition', [0 0 25 25]); 
xlabel('$\mu$','Interpreter','latex')
ylabel('$s_{min}$, RBF','Interpreter','latex')

title(join(["Possible $\mu_T$ for ",curve_param,", $n_{RI}$= ",RI,", N = ",num2str(N)],""),'Interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',23)
filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/transmission/'];
filename = join([filename,curve_param,"signflip",mu_min,"_",mu_max,"_N",num2str(N),"_RBF"],'');
filename = strrep(filename,'.',"dot");
% exportgraphics(gcf,join([filename,".png"],""),'Resolution',300)  
% savefig(gcf,filename)