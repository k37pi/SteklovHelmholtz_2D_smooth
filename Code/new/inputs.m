clear; close all
set(0,'defaultTextInterpreter','latex'); 
set(0,'defaultAxesFontSize',20)
set(0, 'DefaultLineLineWidth', 2);

%% curve names, add curve name if you updated the curve library ---------
curves_list = ["" + ...
    "Circle"; "Ellipse"; "kite";"coolcurve";"disk_def";
    "acorn";"star";"rectangle";"egg";"opt_disk";"butterfly";
    "cavity";"bunimovich";"sponge";"circle2";"sponge2";"coolpolar";
    "pert_circle";"ellipse_opt";"tear";"coolpolar_1d";"polar2";
    "asypolar";"coolpolar_2d";"flower";"insect";"fourier"];%;"peanut"];

% circle = 1, ellipse = 2,  kite = 3, coolcurve = 4, 
% disk_def = 5, acorn = 6 , star = 7, square = 8, egg = 9,
% opt_disk = 10, butterfly = 11, cavity = 12, bunimovich = 13,
% sponge = 14, shifted circle = 15, shifted sponge = 16,
% coolpolar = 17, pert_circle = 18, ellipse_opt = 19, tear = 20;
% coolpolar_1d = 21, polar2 = 22, asypolar = 23, coolpolar_2d = 24;
% flower = 25, insect = 26, fourier = 27
curve_number = 27;

% second ev as fn ecc for fixed mu, many mu

%% mu, N, M ----------
Mu = (7);19.951;4995891589;300.005956478458;  %0.8824
mu = Mu; % wave number for helmholtz
N = 250; % half the number of points
N1 = 45; % half the number of points for init svd of single layer for Alg. 1
M = 150; % number of interior points(ish)

%% curve params --------
disk_r = 1;
ellipse_a = 5;    
ellipse_b = 2; ellipse_ecc = 0.3;
square_a = 1; square_b = 1;

epps_egg = 0; 
aa = -1; bb = 1;
cooln = 7; % have fixed first 9, only last changes.

Minr = zeros(1,2*cooln)-1; % min of each parameter is -0.2
Maxr = zeros(1,2*cooln)+1; % max of each parameter is 0.2

coolp_ab = Minr+(Maxr-Minr).*rand(1,cooln*2);[0.5000   -0.5000    0.1252    0.5000];[0.1 0 0 0 0.1;0 0.1 0 0 -0.1];%coolp_params([end-1 end],:);%;%[-0.1847   -0.1998   -0.1998   -0.1999   -0.1997;-0.1982    0.1996   -0.1997    0.1995   -0.1997];
coolp_a0 = 1;
cool_ab = [Minr+(Maxr-Minr)*rand(2*cooln)];
coolp_ab_al = [-0.48245    0.32807    -0.18202    -0.88027    -0.08484    -0.035442];
coolp_ab_ll = [-0.88608    0.39712    0.50742    -0.9832    0.98525    -0.48008];
% coolp_ab = [0.0015   -0.0003   -0.0015    0.0062   -0.9469    0.0059    0.0026    0.0040    0.0015   -0.0467];
% % coolp_ab = [ -0.0057    0.0056    0.0036   -0.0005    0.0032    0.9588   -0.0072    0.0083    0.0012   -0.0028   -0.0030    0.0010];
%  coolp_ab = [ 0.0005   -0.0008   -0.0013    0.0003   -0.0010    0.0003   -0.0006    1.0000    0.0003   -0.0000    0.0005    0.0003 0.0004   -0.0004   -0.0006    0.0300];
% coolp_ab = [-0.0054   -0.0047   -0.0008   -0.0043    0.0008   -0.0030   -0.9715   -0.0008   -0.0049   -0.0036    0.0054  -0.0155    0.0191   -0.0054];
% coolp_ab=[-0.21813    -0.87769    0.64259    0.41541    -0.39406    0.96825    -0.45356    0.76863    -0.28111    -0.16671];
% coolp_ab = [-3.1111    -4.6193    -1.6395    2.0918    1.1323    2.5676    4.4027    3.9472    -3.1351    -1.3177];

fouriername = "flo5";
asx = readmatrix(join(["files/",fouriername,"_ax.csv"],"")); bsx = readmatrix(join(["files/",fouriername,"_bx.csv"],""));
asy = readmatrix(join(["files/",fouriername,"_ay.csv"],"")); bsy = readmatrix(join(["files/",fouriername,"_by.csv"],""));

diskdef_k = 1.5; 
sa = -0.2; % sponge a

Hom = 1; % domain scaling
len1 = 0; % make perimeter a (make this proper)

%-0.5+1*rand(1)
curve_params = {asx, bsx, asy, bsy};[3 4 2 4];coolp_ab;
% [0.1	-0.1	1.26741544171621E-08	8.25734536375718E-09];[1 10];[0  -0.1   0   0];
% [-0.0000   -0.0000   -0.1000   -0.0000   -0.0000    0.1000];
% [7.45174500679756e-10	-0.1	5.35429151946838e-11	3.16015183222571e-10];
% [-0.0932   -0.0423    0.0993    0.0384];
% 1.0e-08 *[0.1167    0.0169];
% 1e-3*[-0.101313862634610   0.000024669315423  -0.000000940231285   0.000000662229649   0.003070416525896  -0.000001250566647 0.000001150719191  -0.000000475009567];
% [-0.099596996613530   0.029847470617342   0.072455734671286  -0.091657401841007  -0.044123507854692   0.095469917041900];
% [ -0.081939761721186  -0.023237610369665  -0.082659884105312  -0.042107710259504];
% [0.082590127751369  -0.000000000005976];
% 1.0e-03*[-0.3708    0.3708];
% [0.0625726005358623	-0.0201281674546541	0.0951046704728284	-0.0671765026172839	-0.00353359866488979	-0.0999965268476433];
%[0.0626   -0.0201    0.0951   -0.0672   -0.0035   -0.1000];[0.1000    0.0036    0.0690    0.0095];[0.0654   -0.0034    0.0949    0.0090];[0 .1 0 .1];[-0.1453    0.0161    0.9922    0.0317   -0.0843    0.0462];[0.8963   -0.0223   -0.9926   -0.2075];coolp_ab;[ellipse_a ellipse_b];coolp_ab;0.86193; % set this as the parameter for the domain 
%curve_params = 0;coolp_ab(1);%[coolp_ab(1,:) coolp_ab(2,:)]; 

intext = "[]"; % keep intext empty for Steklov-Laplace
%% svd tolerance ------
tol = 1e-4;
svd_flag = 1; % needed only if not using 'alg1_sh.m'

%% eigenvalue numbers ---------
fixes = [1 7]; % [1:2:5 10 15];

%% choice of solution reconstruction ---------
% this is a text string for the 'efn_in' function to select which solution 
% to reconstruct inside the domain. 
% "l" = dirichlet/neumann/steklov-laplace, "h" = dirichlet/neumann/steklov-helmholtz
% "d" = dirichlet-laplace evp
% use this as required in the correct situation, see the section on efn
% plots in 'steklov_HL.m' or 'helm_bvp.m'
LH = "l"; 

%% shape opt --------
k = 1; mult = 1; pm = 1;

%% points and all that [ no input required =) ] ------
[t,x,dx,d2x,nx,points_inside,meshes,X,R,dp1norm,dp2norm,len,dom_area,curve_name] = point_setup(curve_number,curve_params,N,M,len1,Hom);
kappa = (dx(1,1:end-1).*d2x(2,1:end-1)-d2x(1,1:end-1).*dx(2,1:end-1))./(dx(1,1:end-1).^2+dx(2,1:end-1).^2).^1.5;
kappa = sum(kappa)*pi/N; kappa = kappa/(2*pi);
%
figure(1)
plot(x(1,:),x(2,:)); hold on
quiver(x(1,1),x(2,1),nx(1,1),nx(2,1),'LineWidth',2.5,'MaxHeadSize',0.5,'color',"#77AC30")
axis equal

%% boundary curve_paramsweights --------
rho_wt = @(t) 1;exp(cos(t));
%rho_wt = @(x,nx) 1./(sum(x.*nx));
if nargin(rho_wt)==1
Rho_wt = rho_wt(t);
else
Rho_wt = rho_wt(x,nx);
end

% figure;
% plot3(x(1,:),x(2,:),Rho_wt,'--*')
%figure
% colormap(flipud(gray))