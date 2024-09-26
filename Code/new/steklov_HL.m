% -----------------------------------------------------------------------% 
% This script computes the Steklov spectrum for the 
% Laplace and Helmholtz equations in the plane for smooth enough domains. 
% The inputs are set in the 'inputs.m' file (duh doi). 

% These scripts are required:
% ftp.m || curve.m || point_setup.m || stek_helm.m || alg1_sh.m || stek_lap.m
% efn_in.m || inputs.m 

% The sections in this file are:
% *load parameters - it loads all the good stuff you set from inputs.m
% *solvers - the alg1 and stek_lap functions return the 
%  evs, density, functions and layer pot matrices.
% *the other 2 sections are plots for evs and efns. 

% you can use the spectrum to do whatever you want with it. 
% -----------------------------------------------------------------------% 
clear; close all
set(0,'defaultTextInterpreter','latex'); set(0,'defaultAxesFontSize',20)

%% load parameters -------
run 'inputs.m'

scale_mu = mu;%/sqrt(dom_area); 
evk = 0;

%% solvers --------
[rank_def,evs_sh,eden_sh,evecs_sh,A_sh,B_sh] = alg1_sh(curve_number,curve_params,scale_mu,N1,N,M,len1,tol,Hom,rho_wt,evk);
[evs_sl,eden_sl,evecs_sl,A_sl,B_sl] = stek_lap(curve_number,curve_params,N,M,len1,Hom,rho_wt,intext);

%% eigenvalue plots ------ 
k = 1:N; % which evs you want to see ?
figure
plot(k,evs_sl(k),'--*')
hold on
plot(k,evs_sh(k),'--*')
xlabel("$k$"); ylabel("$\sigma_k$")
legend(["SL","SH"],'Interpreter','latex')

cn = strrep(curve_name,"_","\_");

if curve_name ~= "fourier"
title_text = join([cn,curve_params(1),", $\mu=$ ",num2str(Mu),", rank deficiency SH = ",rank_def],"");
else
title_text = join([cn,", $\mu=$ ",num2str(Mu),", rank deficiency SH = ",rank_def],"");
end
title(title_text)

figure()
plot(evs_sh(1:N)-evs_sl(1:N),'--*')
xlabel("$k$"); ylabel("$\sigma^{SL}_k-\sigma^{SH}_k$")
title(join(["$\mu$ = ",mu],""))

%% eigenfunction plots ------ 
% fixes is set in 'inputs.m' but, you can edit it here too. 
% For ex here its the first ev, the first positive SH ev and the 25th.
% fixes = [1, find(evs_sh > 0,1), 93]; 
fixes = [1:5]; 

[evec_in_sl,Xsl,Ysl] = efn_in(fixes,eden_sl,x,dx,points_inside,0,M,"l");
[evec_in_sh,Xsh,Ysh] = efn_in(fixes,eden_sh,x,dx,points_inside,Mu,M,"h");

for i = 1:length(fixes)    
figure(10)
ax1 = subplot(2,length(fixes),i);
surf(Xsl{i},Ysl{i},real(evec_in_sl{i}));
shading interp; view(0,90); colormap(ax1,viridis); %colorbar
title(join(["$\sigma^{SL}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sl(fixes(i))) ],"") )
xticks(nan); yticks(nan) 
axis image
axis tight

ax2 = subplot(2,length(fixes),i+length(fixes));
surf(Xsh{i},Ysh{i},real(evec_in_sh{i}));
shading interp; view(0,90); colormap(ax2,magma); %colorbar
title(join(["$\sigma^{SH}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sh(fixes(i))) ],"") )
xticks(nan); yticks(nan) 

axis image
axis tight

figure(11)
ax1 = subplot(2,length(fixes),i);
surf(Xsl{i},Ysl{i},abs(evec_in_sl{i}).^2);
shading interp; view(0,90); colormap(ax1,viridis); %colorbar
title(join(["$\sigma^{SL}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sl(fixes(i))) ],"") )
xticks(nan); yticks(nan) 
axis image
axis tight

ax2 = subplot(2,length(fixes),i+length(fixes));
surf(Xsh{i},Ysh{i},abs(evec_in_sh{i}).^2);
shading interp; view(0,90); colormap(ax2,magma); %colorbar
title(join(["$\sigma^{SH}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sh(fixes(i))) ],"") )
xticks(nan); yticks(nan) 

axis image
axis tight

end
sgt_text = join([cn,", $|M|$ = ",len,", $\mu=$ ",num2str(Mu)],"");
if rank_def ~=0 
   sgt_text = join([sgt_text,", SVD"],""); 
end    
figure(10)
sgtitle(join([sgt_text,", $Re(u)$"],""),'fontsize',23)

figure(11)
sgtitle(join([sgt_text,", $|u|^2$"],""),'fontsize',23)

%% potential matrices 
figure()

MAR = max(max(real(A_sh))); mAR = min(min(real(A_sh)));
MBR = max(max(real(B_sh))); mBR = min(min(real(B_sh)));

MAI = max(max(imag(A_sh))); mAI = min(min(imag(A_sh)));
MBI = max(max(imag(B_sh))); mBI = min(min(imag(B_sh)));

subplot(2,4,1)
imagesc(real(A_sh)); %colorbar
title("real A")
axis image
clim([mAR MAR])

subplot(2,4,2)
imagesc(imag(A_sh)); %colorbar
title("imag A")
axis image
clim([mAI MAI])

subplot(2,4,3)
imagesc(real(A_sh-A_sh')); %colorbar
title("real(A-A')")
axis image
clim([mAR MAR])

subplot(2,4,4)
imagesc(imag(A_sh-A_sh')); %colorbar
title("imag(A-A')")
axis image
clim([mAI MAI])
colorbar

subplot(2,4,5)
imagesc(real(B_sh)); %colorbar
title("real B")
axis image
clim([mBR MBR])

subplot(2,4,6)
imagesc(imag(B_sh)); %colorbar
title("imag B")
axis image
clim([mBI MBI])

subplot(2,4,7)
imagesc(real(B_sh-B_sh')); %colorbar
title("real(B-B')")
axis image
clim([mBR MBR])

subplot(2,4,8)
imagesc(imag(B_sh-B_sh')); 
title("imag(B-B')")
axis image
clim([mBI MBI])
colorbar

%% repeated eigs ------------------
[~,blah,ev_reps] = uniquetol(real(evs_sh));
ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals