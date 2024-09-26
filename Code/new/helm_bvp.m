% -----------------------------------------------------------------------% 
% This script solves Dirichlet, Neumann BVPs for
% Laplace and Helmholtz equations in the plane for smooth enough domains. 
% Note: Neumann-Laplace is better conditioned for exterior problem.  
% The inputs are set in the 'inputs.m' file (duh doi). 

% These scripts are required:
% ftp.m || curve.m || point_setup.m || stek_helm.m || stek_lap.m
% efn_in.m || inputs.m 

% -----------------------------------------------------------------------% 
clear; close all
set(0,'defaultTextInterpreter','latex'); set(0,'defaultAxesFontSize',20)
    

%% load basic parameters from inputs.m ------------------
run inputs.m 

%% extra inputs: select boundary term, other wave number ----------------
f = @(t) sin(2*t);%-cos(t);%(-1).^floor(t); 
%f = @(x) x(1,:).^2.*x(2,:);
int1 = integral(f,0,2*pi);
if abs(int1) > 1e-6
disp("f does not integrate to 0 on the boundary !")
end
f_latex = "$\sin 2t$";%"$-1^{|\_t\_|}$";
f_grid = transpose(f(t(1:end-1)));
%f_grid = transpose(f(x(:,1:end-1)));

% uex = @(x) ([1 0]*x) .* ([0 1]*x-0.3); 
% f_grid = uex(x); f_grid = transpose(f_grid(1:end-1));

mu1 = 0.7;
%% BVPs --------
rho_wt = @(t) 1; % keep this as 1 for this problem
[~,~,~,~,Ah,Bh] = alg1_sh(curve_number,curve_params,mu,N1,N,M,len1,tol,Hom,rho_wt,0);
hdir = Bh\(2*f_grid);
hneum = Ah\(2*f_grid);

[~,~,~,Al,Bl] = stek_lap(curve_number,curve_params,N,M,len1,Hom,rho_wt,[]);
[u,s,v] = svd(Al); ss = diag(s);
ldir = Bl\(2*f_grid);
lneum = Al\(2*f_grid);
P = diag(Al);
%P(P < 1) = 1./P(P<1);
% P(P > 1) = 1./P(P>1);
% lneum = (diag(P)/Al)/(2*f_grid);

[~,x,dx,~,~,points_inside,~,~,~,~,~,len,dom_area,~] = point_setup(curve_number,curve_params,N,M,len1,Hom);
[hdir_in,Xdh,Ydh] = efn_in(1,hdir,x,dx,points_inside,mu,M,"h");
[hneum_in,Xnh,Ynh] = efn_in(1,hneum,x,dx,points_inside,mu,M,"h");
[ldir_in,Xdl,Ydl] = efn_in(1,ldir,x,dx,points_inside,mu,M,"l");
[lneum_in,Xnl,Ynl] = efn_in(1,lneum,x,dx,points_inside,mu,M,"l");

ax1 = subplot(2,3,1);
surf(Xdl{1},Ydl{1},real(ldir_in{1}));
shading interp; view(0,90); colormap(ax1,magma); colorbar
hold on 
plot3(x(1,1:end-1),x(2,1:end-1),f_grid,'k')
title("Dirichlet Laplace") 
axis equal

ax2 = subplot(2,3,2);
surf(Xdh{1},Ydh{1},real(hdir_in{1}));
shading interp; view(0,90); colormap(ax2,magma); colorbar
hold on 
plot3(x(1,1:end-1),x(2,1:end-1),f_grid,'k')
title(join(["Dirichlet Helmholtz, $\mu$ = ",mu],"")) 
axis equal

ax4 = subplot(2,3,4);
surf(Xnl{1},Ynl{1},real(lneum_in{1}));
shading interp; view(0,90); colormap(ax4,viridis); colorbar
hold on 
plot3(x(1,1:end-1),x(2,1:end-1),f_grid,'k')
title("Neumann Laplace") 
axis equal


ax5 = subplot(2,3,5);
surf(Xnh{1},Ynh{1},real(hneum_in{1}));
shading interp; view(0,90); colormap(ax5,viridis); colorbar
hold on 
plot3(x(1,1:end-1),x(2,1:end-1),f_grid,'k')
title(join(["Neumann Helmholtz, $\mu$ = ",mu],"")) 
axis equal

[~,~,~,~,Ah,Bh] = alg1_sh(curve_number,curve_params,mu1,N1,N,M,len1,tol,Hom,rho_wt,0);

hdir = Bh\(2*f_grid);
hneum = Ah\(2*f_grid);
[hdir_in,Xdh,Ydh] = efn_in(1,hdir,x,dx,points_inside,mu1,M,"h");
[hneum_in,Xnh,Ynh] = efn_in(1,hneum,x,dx,points_inside,mu1,M,"h");

ax3 = subplot(2,3,3);
surf(Xdh{1},Ydh{1},real(hdir_in{1}));
shading interp; view(0,90); colormap(ax3,magma); colorbar
hold on 
plot3(x(1,1:end-1),x(2,1:end-1),f_grid,'k')
title(join(["Dirichlet Helmholtz, $\mu$ = ",mu1],"")) 
axis equal

ax6 = subplot(2,3,6);
surf(Xnh{1},Ynh{1},real(hneum_in{1}));
shading interp; view(0,90); colormap(ax6,viridis); colorbar
hold on 
plot3(x(1,1:end-1),x(2,1:end-1),2*f_grid,'k')
title(join(["Neumann Helmholtz, $\mu$ = ",mu1],"")) 
axis equal
%sgtitle(join(["$f(t)$ = ",f_latex, ", $\mu$ = ",mu ],""),'fontsize',25,'interpreter','latex')
sgtitle(join(["$f(t)$ = ",f_latex],""),'fontsize',25,'interpreter','latex')