% This function returns domain details and required matrices to approx. the
% integrals. 

function [t,x,dx,d2x,nx,points_inside,meshes,X,R,dp1norm,dp2norm,len,dom_area,curve_name] = point_setup(curve_number,curve_params,N,M,len1,Hom)

curves_list = ["" + ...
    "Circle"; "Ellipse"; "kite";"coolcurve";"disk_def";
    "acorn";"star";"rectangle";"egg";"opt_disk";"butterfly";
    "cavity";"bunimovich";"sponge";"circle2";"sponge2";"coolpolar";
    "pert_circle";"ellipse_opt";"tear";"coolpolar_1d";"polar2";
    "asypolar";"coolpolar_2d";"flower";"insect";"fourier"];
% cool curve = (a1 cos (b1 t) + a2 sin (b2 t))ang2a[cost, sint] 
% circle = 1, ellipse = 2,  kite = 3, coolcurve = 4, 
% disk_def = 5, acorn = 6 , star = 7, square = 8, egg = 9,
% opt_disk = 10, butterfly = 11, cavity = 12, bunimovich = 13,
% sponge = 14, shifted circle = 15, shifted sponge = 16,
% coolpolar = 17, pert_circle = 18, ellipse_opt = 19, tear = 20
% coolpolar_1d = 21, polar2 = 22, asypolar = 23, coolpolar_2d = 24
% flower = 25
curve_name = curves_list(curve_number); curve_name = lower(curve_name);

% curve params
hom = Hom; % domain scaling

% inside grid
m = M;
%--------------------------------------------------------

%% curve points and parametrization ---------------------
a = 0; b = 2*pi;
dt = (b-a)/(2*N); %cutoff = 0.1; 
t = a:dt:b; % 2N+1 points

switch curve_name
    case "circle"
        disk_r = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,disk_r);
        V = [cos(t(1:end-1));sin(t(1:end-1))]; V = {V};
        %curve_param = join([curve_name,disk_r],"");
    case "ellipse"
        ellipse_a = curve_params(1); 
        ellipse_b = curve_params(2);

        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,ellipse_a,ellipse_b);
        %curve_param = join([curve_name,num2str(ellipse_a),'\_',num2str(ellipse_b)],"");
    case "ellipse_opt"
        ellipse_a = curve_params(1); 

        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,ellipse_a);
        %curve_param = join([curve_name,num2str(ellipse_a),'\_',num2str(ellipse_b)],"");
  
    case "coolcurve"
        cool_ab = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,cool_ab);
        while(~isempty(selfintersect(x(1,1:end-1),x(2,1:end-1))))
            break
            % aa = -1/4; bb = -aa;
            % cool_ab = aa+(bb-aa)*rand(4,1);
            % [x,dx,d2x,nx,len] = curve(curve_name,t,len1,cool_ab);
        end    
    case "coolpolar"
        coolp_ab = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,coolp_ab,1);
        %curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3)),'\_',num2str(cool_ab(4)),'\_'],"");
        %curve_param = join([curve_name,num2str(coolp_ab(1))],"");
        % while(~isempty(selfintersect(x(1,1:end-1),x(2,1:end-1))))
        %     coolp_ab = coolp_ab+0.01;
        %     [x,dx,d2x,nx,len] = curve(curve_name,t,len1,coolp_ab,1);
        %     %curve_param = join([curve_name,num2str(coolp_ab(1)),'\_',num2str(coolp_ab(2)),'\_',num2str(coolp_ab(3))],"");
        % end      
    case "coolpolar_1d"
        coolp_ab = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,coolp_ab,1);
        %curve_param = join([curve_name,num2str(cool_ab(1)),'\_',num2str(cool_ab(2)),'\_',num2str(cool_ab(3)),'\_',num2str(cool_ab(4)),'\_'],"");
        %curve_param = join([curve_name,num2str(coolp_ab(1))],"");
        % while(~isempty(selfintersect(x(1,1:end-1),x(2,1:end-1))))
        %     coolp_ab = coolp_ab+0.01;
        %     [x,dx,d2x,nx,len] = curve(curve_name,t,len1,coolp_ab,1);
        %     %curve_param = join([curve_name,num2str(coolp_ab(1)),'\_',num2str(coolp_ab(2)),'\_',num2str(coolp_ab(3))],"");
        % end      
    case "coolpolar_2d"
        coolp_ab = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,coolp_ab,1);
    case "disk_def"
        diskdef_k = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,diskdef_k);  
        %curve_param = join([curve_name,diskdef_k],"");
    case "rectangle"
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,square_a,square_b);  
        %curve_param = join([curve_name,num2str(square_a),'\_',num2str(square_b)],"");
      case "egg"
        epps_egg = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,epps_egg);  
        %curve_param = join([curve_name,num2str(epps_egg)],"");
    case "sponge"
        sa = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,sa);  
    case "flower"
        fl = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,fl);  
    case "insect"
        ins = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,ins);  
    case "fourier"
        asx = curve_params{1}; bsx = curve_params{2};
        asy = curve_params{3}; bsy = curve_params{4};
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,asx,bsx,asy,bsy);  
    case "polar2" 
        a2 = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,a2);  
    case "asypolar"
        a2 = curve_params;
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1,a2);  
    otherwise
        [x,dx,d2x,nx,len] = curve(curve_name,t,len1);
        %curve_param = curve_name;
end 


if curve_name ~= "circle" && len1 ~= 1
    x = hom*x; dx = hom*dx; d2x = hom*d2x; len = hom*len;
    nx = [dx(2,:);-dx(1,:)];
    nx = nx./sqrt(dx(1,:).^2+dx(2,:).^2);
end

dom_area = 2*pi/(length(t)-1)*sum((dx(1,1:end-1).*x(2,1:end-1)));%polyarea(x(1,:)',x(2,:)');
dom_area = abs(dom_area);
kappa = (dx(1,1:end-1).*d2x(2,1:end-1)-d2x(1,1:end-1).*dx(2,1:end-1))./(dx(1,1:end-1).^2+dx(2,1:end-1).^2).^1.5;
kappa = sum(kappa)*pi/N; kappa = kappa/(2*pi);

xi = linspace(min(x(1,:)),max(x(1,:)),m);
yi = linspace(min(x(2,:)),max(x(2,:)),m);
[X,Y] = meshgrid(xi,yi);
idx = inpolygon(X(:),Y(:),x(1,:),x(2,:)) ;
points_inside = transpose([X(idx),Y(idx)]);

points_inside_close = 0.99*x(:,1:end-1);

points_inside = [points_inside,points_inside_close];

[Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
[Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

[Tau,T] = meshgrid(t);
[Xtau,Xt] = meshgrid(x(1,:));
[Ytau,Yt] = meshgrid(x(2,:));
[dXtau,dXt] = meshgrid(dx(1,:));
[dYtau,dYt] = meshgrid(dx(2,:)); 
[d2Xtau,d2Xt] = meshgrid(d2x(1,:));
[d2Ytau,d2Yt] = meshgrid(d2x(2,:));

meshes = {Xtauin,Xin,Ytauin,Yin,Tau,T,Xtau,Xt,Ytau,Yt,dXtau,dXt,dYtau,dYt,d2Xtau,d2Xt,d2Ytau,d2Yt};

X = dYt.*(Xtau-Xt)-dXt.*(Ytau-Yt);
R = sqrt((Xt-Xtau).^2+(Yt-Ytau).^2);

dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);

end