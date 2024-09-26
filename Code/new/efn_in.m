function [evec_in,X,Y] = efn_in(fixes,eden,x,dx,points_inside,mu,M,LH)
m = M;
evec_in = cell(length(fixes),1);
X = cell(length(fixes),1);
Y = cell(length(fixes),1);
N = size(eden,1)/2;
t = 0:pi/N:(2*pi); [~,T] = meshgrid(t);
t_scaled = (0:2*N).*T; % the first column is 0*t, second is 1*t, so on.
cos_vecs = cos(t_scaled); % l is fixed in columns, t_j is fixed in rows.
sin_vecs = sin(t_scaled);

xi = linspace(min(x(1,:)),max(x(1,:)),M);
yi = linspace(min(x(2,:)),max(x(2,:)),M);
[Xs,Ys] = meshgrid(xi,yi);
points_inside = transpose([Xs(:),Ys(:)]);

points_inside_close = 0.99*x(:,1:end-1);
points_inside = [points_inside,points_inside_close];

[Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
[Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
Rin = sqrt(Xindiff.^2+Yindiff.^2);

[dXtau,~] = meshgrid(dx(1,:));
[dYtau,~] = meshgrid(dx(2,:)); 

dp2norm = sqrt(dXtau.^2+dYtau.^2);

for ff = 1:length(fixes)
fix = fixes(ff);
eden_fix = eden(:,fix);

if lower(LH) == "h"
HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
Recon_integrand = HankRin.*transpose(eden_fix).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); 
elseif lower(LH) == "l"
HankRin = log(Rin); % rows correspond to fixed point inside 
Recon_integrand = HankRin.*transpose(eden_fix).*dp2norm(fix,(1:end-1)); % make this a tp and then sum and scale
recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
recint_sin = recint_sin1*sin_vecs(:,2:N)';
recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
solution_fixed_point = -1/(2*N)*sum(recint,2); 
elseif lower(LH) == "d" || lower(LH) == "n"
    [~,S,V] = svd(eden); %toc
    S = diag(S); [~,min_ind] = min(S); 
    %all_S(:,ms) = S;
    phi_min = V(:,min_ind); %evec_dir = B/2*phi_min;
    HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
    Recon_integrand = HankRin.*transpose(phi_min).*dp2norm(min_ind,(1:end-1)); % make this a tp and then sum and scale
    recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
    recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2; 
    recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
    recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
    recint_sin = recint_sin1*sin_vecs(:,2:N)';
    recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
    solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
% else
%     HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
%     Recon_integrand = HankRin.*transpose(eden).*dp2norm(1,(1:end-1)); % make this a tp and then sum and scale
%     recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
%     recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2; 
%     recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
%     recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
%     recint_sin = recint_sin1*sin_vecs(:,2:N)';
%     recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
%     solution_fixed_point = 1i*pi/(4*N)*sum(recint,2); % trapz rule, by matlab below
end

xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));

[Xin2,Yin2] = meshgrid(xlin,ylin);

evec_in{ff} = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;
X{ff} = Xin2; Y{ff} = Yin2;
end

end