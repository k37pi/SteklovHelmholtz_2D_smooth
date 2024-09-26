function [evs,eden,evecs,A,B] = stek_lap(curve_number,curve_params,N,M,len1,Hom,rho_wt,intext)
[t,x,~,~,nx,~,meshes,X,R,dp1norm,dp2norm,~,~,curve_name] = point_setup(curve_number,curve_params,N,M,len1,Hom);

Tau = meshes{5}; T = meshes{6};
dXt = meshes{12};
dYt = meshes{14};
d2Xt = meshes{16};
d2Yt = meshes{18};

%% R, X for all pairs t, tau ------
XR2_ratio = X./(R.^2);

%% L, M ---------------------------------------
dpnorm_ratio = dp2norm./dp1norm;

logsin_term = log(4*sin((T-Tau)/2).^2);
cL = (1/pi)*dpnorm_ratio;
L = cL.*XR2_ratio; 
L(end,1) = 0; L(1,end) = 0;

if curve_name == "rectangle"
    L(isnan(L)) = 0; 
else
    L(isnan(L)) = diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2)); 
end
L2 = L;
    
% there is a negative sign here, but we just multiply the 
% evs with -1 at the end
M = (1/pi)*dp2norm.*log(R);
M(1,end) = 0; M(end,1) = 0;
M1 = 1/(2*pi)*dp2norm;

M2 = M-M1.*logsin_term;
M2(1,end) = 0; M2(end,1) = 0;

if curve_name == "rectangle"
    M2(isnan(M2)) = 0;
else    
    %M2(isnan(M2)) = diag(dp1norm.*(mc-log(mu/2*dp1norm)/pi));
    M2(isnan(M2)) = diag(dp1norm.*log(dp1norm)/pi);
end

%% Trig. poly. -----------------------------
% f^tp(t,tau) = 
% (a_0/2,a_1,...,a_{N-1},a_N/2).(1,cos(tau),...,cos([N-1]tau),cos([N]tau))+
% (b_1,...,b_{N-1}).(sin(tau),...,sin([N-1]tau)) =
% a_vec.cos_0N_vec + b_vec.sin_0N_vec
t_scaled = (0:2*N).*T; % the first column is 0*t, second is 1*t, so on.
cos_vecs = cos(t_scaled); % l is fixed in columns, t_j is fixed in rows.
sin_vecs = sin(t_scaled);

%L1 tp
%L1_tp = ftp(L1,N,cos_vecs,sin_vecs);

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

RN_mat = cos(diff_t_tau);
for k1 = 2:(N-1)
    RN_mat = RN_mat+cos(k1*diff_t_tau)/k1; 
end    
RN_mat = -2*pi*RN_mat/N - pi*cos(N*(diff_t_tau))/N^2;

%% Matrices ------------------------------
A = pi/N*L2_tp; % (A + I) is the effect of K'
A(:,end) = []; A(end,:) = [];

I = diag(ones(2*N,1));

B = RN_mat.*M1_tp+pi/N*M2_tp; % B is the effect of the S 
B(:,end) = []; B(end,:) = [];
B = -B;

if nargin(rho_wt)==1
Rho_wt = rho_wt(t(1:end-1))';
else
Rho_wt = rho_wt(x(:,1:end-1),nx(:,1:end-1))';
end
B_wt = Rho_wt.*B;

%% Steklov eigenstuff --------------------------------
[eden1,evs1] = eig(I+A,B_wt); 

B = B_wt; 
if isempty(intext)
A = A+I;
else
A = A-I;
end

evs1 = diag(evs1);
evs = real(evs1);    
[evs, ind] = sort(evs);%,"descend");
eden = eden1(:,ind); % eigen densities
evecs = B_wt*eden/2; %disk_r

end