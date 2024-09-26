function [rank_def,evs,eden,evecs,A,B] = stek_helm(curve_number,curve_params,Mu,N,M,len1,tol,svd_flag,Hom,rho_wt,eigsk)
mu = Mu;
[t,x,~,~,nx,~,meshes,X,R,dp1norm,dp2norm,~,~,curve_name] = point_setup(curve_number,curve_params,N,M,len1,Hom);

% Xtauin = meshes{1}; Xin = meshes{2};
% Ytauin = meshes{3}; Yin = meshes{4};
Tau = meshes{5}; T = meshes{6};
% Xtau = meshes{7}; Xt = meshes{8};
% Ytau = meshes{9}; Yt = meshes{10};
%dXtau = meshes{11}; 
dXt = meshes{12};
%dYtau = meshes{13}; 
dYt = meshes{14};
%d2Xtau = meshes{15}; 
d2Xt = meshes{16};
%d2Ytau = meshes{17}; 
d2Yt = meshes{18};

%% R, X for all pairs t, tau
XR_ratio = X./R;

%% L, M ---------------------------------------
dpnorm_ratio = dp2norm./dp1norm;

logsin_term = log(4*sin((T-Tau)/2).^2);
cL = (1i*mu/2)*dpnorm_ratio;
cL1 = -mu/(2*pi)*dpnorm_ratio;
L = cL.*XR_ratio.*besselh(1,mu*R); % note that besselh(1,0) blows up so when R = 0 it's useless.
L(end,1) = 0; L(1,end) = 0;

L1 = cL1.*XR_ratio.*besselj(1,mu*R);

L2 = L-L1.*logsin_term;
L2(end,1) = 0; L2(1,end) = 0;

if curve_name == "rectangle"
    %L(isnan(L)) = 0; 
    L2(isnan(L2)) = 0; 
else
    L(isnan(L)) = diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2)); 
    L2(isnan(L2)) = diag(L); 
end
    
L1(isnan(L1)) = 0;


M = (1i/2)*dp2norm.*besselh(0,mu*R);
M(1,end) = 0; M(end,1) = 0;
M1 = -1/(2*pi)*dp2norm.*besselj(0,mu*R);

M2 = M-M1.*logsin_term;
M2(1,end) = 0; M2(end,1) = 0;
mc = 1i/2-double(eulergamma)/pi;

% if Mu > 0
if curve_name == "rectangle"
    M2(isnan(M2)) = 0;
else    
    M2(isnan(M2)) = diag(dp1norm.*(mc-log(mu/2*dp1norm)/pi));
end

t_scaled = (0:2*N).*T; % the first column is 0*t, second is 1*t, so on.
cos_vecs = cos(t_scaled); % l is fixed in columns, t_j is fixed in rows.
sin_vecs = sin(t_scaled);

%L1 tp
L1_tp = ftp(L1,N,cos_vecs,sin_vecs);

% L2 tp
L2_tp = ftp(L2,N,cos_vecs,sin_vecs);

%M1 tp
M1_tp = ftp(M1,N,cos_vecs,sin_vecs);

% M2 tp
M2_tp = ftp(M2,N,cos_vecs,sin_vecs);

%% RN weights, same for any integrand (L1,M1 for us) ---------------
tic
diff_t_tau = T-Tau;

RN_mat = cos(diff_t_tau);
for k1 = 2:(N-1)
    RN_mat = RN_mat+cos(k1*diff_t_tau)/k1; 
end    
RN_mat = -2*pi*RN_mat/N - pi*cos(N*(diff_t_tau))/N^2;

%% Matrices ------------------------------
A = RN_mat.*L1_tp+pi/N*L2_tp; % (A + I) is the effect of K'
A(:,end) = []; A(end,:) = [];

I = diag(ones(2*N,1));

B = RN_mat.*M1_tp+pi/N*M2_tp; % B is the effect of the S 
B(:,end) = []; B(end,:) = [];

if nargin(rho_wt)==1
Rho_wt = rho_wt(t(1:end-1))';
else
Rho_wt = rho_wt(x(:,1:end-1),nx(:,1:end-1))';
end
B_wt = Rho_wt.*B;

%% Steklov eigenstuff --------------------------------
if svd_flag 
[u,s,v] = svd(B_wt);
ss = diag(s);
rank_def = sum(abs(ss-tol) < tol);
A_svd = u(:,1:end-rank_def)'*(A+I)*v(:,1:end-rank_def); 
B_wt1 = u(:,1:end-rank_def)'*B_wt*v(:,1:end-rank_def);    
if eigsk == 0
%disp("eig")    
[eden,evs] = eig(A_svd,B_wt1);    
eden = v(:,1:end-rank_def)*eden;
[evs, ind] = sort(real(diag(evs)));
eden = eden(:,ind); 
evecs = B_wt*eden/2;
else
%disp("eigs")    
[eden,evs] = eigs(A_svd,B_wt1,eigsk,"smallestreal");    
eden = v(:,1:end-rank_def)*eden;
[evs, ind] = sort(real((evs)));
eden = eden(:,ind); 
evecs = B_wt*eden/2;    
end    
A = A_svd; B = B_wt1;
else
if eigsk == 0    
%disp("eig")    
[eden,evs] = eig(I+A,B_wt); 
[evs, ind] = sort(real(diag(evs)));
eden = eden(:,ind);
evecs = B_wt*eden/2;
else
%disp("eigs")    
[eden,evs] = eigs(I+A,B_wt,eigsk,"smallestreal"); 
[evs, ind] = sort(real(diag(evs)));
eden = eden(:,ind);
evecs = B_wt*eden/2;
end    
A = I+A; B = B_wt;
rank_def = 0;
end

end