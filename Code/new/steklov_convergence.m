% -----------------------------------------------------------------------% 
% This script computes the convergence of the Steklov spectrum 
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
set(0,'DefaultLineLineWidth',2)

%% load parameters -------
run 'inputs.m'

%% extra inputs --------
N2 = 1000; % half the number of points for 'true' solution 
Ns = 2.^(0:4)*50; % 1/2 # of points to test at 
dts = pi./Ns;

%% evs matrix ---------
Evs_sl = zeros(N2,length(Ns)+1); Evs_sh = zeros(N2,length(Ns)+1);
AEs_sl = zeros(N2,length(Ns)); REs_sl = zeros(N2,length(Ns));
AEs_sh = zeros(N2,length(Ns)); REs_sh = zeros(N2,length(Ns));
ConvRtAE_sl = zeros(N2,length(Ns)-1); ConvRtAE_sh = zeros(N2,length(Ns)-1);
ConvRtRE_sl = zeros(N2,length(Ns)-1); ConvRtRE_sh = zeros(N2,length(Ns)-1);

%% errors --------
disp("Computing true for SH, this will take some time")
[~,evs_sh_true] = alg1_sh(curve_number,curve_params,Mu,N1,N2,M,len1,tol,Hom,rho_wt,0);
disp("Computing true for SL, this will take some time")
[evs_sl_true,~,~,~,~] = stek_lap(curve_number,curve_params,N2,M,len1,Hom,rho_wt,intext);
disp("Longest part done")
evs_sl_true(evs_sl_true < 0) = [];
evs_sl_true(abs(evs_sl_true) < 1e-8) = [];
Evs_sl(1:N2,end) = evs_sl_true(1:N2);
Evs_sh(1:N2,end) = evs_sh_true(1:N2);

for i = 1:length(Ns)
[rank_def,evs_sh,eden_sh,evecs_sh,A_sh,B_sh] = alg1_sh(curve_number,curve_params,Mu,N1,Ns(i),M,len1,tol,Hom,rho_wt,0);
[evs_sl,eden_sl,evecs_sl,A_sl,B_sl] = stek_lap(curve_number,curve_params,Ns(i),M,len1,Hom,rho_wt,intext);
evs_sl(evs_sl < 0) = [];
evs_sl(abs(evs_sl) < 1e-8) = [];
Evs_sl(1:Ns(i),i) = evs_sl(1:Ns(i));
Evs_sh(1:Ns(i),i) = evs_sh(1:Ns(i));
AEs_sl(1:Ns(i),i) = abs(Evs_sl(1:Ns(i),i)-Evs_sl(1:Ns(i),end));
REs_sl(1:Ns(i),i) = AEs_sl(1:Ns(i),i)./abs(Evs_sl(1:Ns(i),end));
AEs_sh(1:Ns(i),i) = abs(Evs_sh(1:Ns(i),i)-Evs_sh(1:Ns(i),end));
REs_sh(1:Ns(i),i) = AEs_sh(1:Ns(i),i)./abs(Evs_sh(1:Ns(i),end));
if i > 1
ConvRtAE_sl(1:Ns(i-1),i-1) = abs(log(AEs_sl(1:Ns(i-1),i)./AEs_sl(1:Ns(i-1),i-1))./log(dts(i)/dts(i-1)));
ConvRtRE_sl(1:Ns(i-1),i-1) = abs(log(REs_sl(1:Ns(i-1),i)./REs_sl(1:Ns(i-1),i-1))./log(dts(i)/dts(i-1)));
ConvRtAE_sh(1:Ns(i-1),i-1) = abs(log(AEs_sh(1:Ns(i-1),i)./AEs_sh(1:Ns(i-1),i-1))./log(dts(i)/dts(i-1)));
ConvRtRE_sh(1:Ns(i-1),i-1) = abs(log(REs_sh(1:Ns(i-1),i)./REs_sh(1:Ns(i-1),i-1))./log(dts(i)/dts(i-1)));
figure(2)
subplot(1,2,1)
semilogy(1:Ns(i-1),ConvRtAE_sl(1:Ns(i-1),i-1),'--*')
xlabel("$k$"); ylabel("$Rate AE_k$")
title("Laplace")
hold on

% subplot(2,2,2)
% semilogy(1:Ns(i-1),ConvRtRE_sl(1:Ns(i-1),i-1),'--*')
% xlabel("$k$"); ylabel("$RE_k$")
% hold on
% title("Laplace")

subplot(1,2,2)
semilogy(1:Ns(i-1),ConvRtAE_sh(1:Ns(i-1),i-1),'--*')
xlabel("$k$"); ylabel("$Rate AE_k$")
title("Helmholtz")
hold on

% subplot(2,2,4)
% semilogy(1:Ns(i-1),ConvRtRE_sh(1:Ns(i-1),i-1),'--*')
% xlabel("$k$"); ylabel("$RE_k$")
% hold on
% title("Helmholtz")
end

figure(1)
subplot(2,2,1)
semilogy(1:Ns(i),AEs_sl(1:Ns(i),i),'--*')
xlabel("$k$"); ylabel("$AE_k$")
title("Laplace")
hold on

subplot(2,2,2)
semilogy(1:Ns(i),REs_sl(1:Ns(i),i),'--*')
xlabel("$k$"); ylabel("$RE_k$")
hold on
title("Laplace")

subplot(2,2,3)
semilogy(1:Ns(i),AEs_sh(1:Ns(i),i),'--*')
xlabel("$k$"); ylabel("$AE_k$")
title("Helmholtz")
hold on

subplot(2,2,4)
semilogy(1:Ns(i),REs_sh(1:Ns(i),i),'--*')
xlabel("$k$"); ylabel("$RE_k$")
hold on
title("Helmholtz")

disp(round(i/length(Ns),2)*100)
end

figure(1)
leg = legend(string(2*Ns),'Interpreter','latex','linewidth',2);
title(leg,"$N$")
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.92;
Lgnd.Position(2) = 0.4;

sgtitle_text = join([curve_name,", $|M|$=",len,", $\mu=$ ",num2str(Mu),", 'true' $N$=",2*N2],"");
if rank_def ~=0
sgtitle_text = join([sgtitle_text,", SVD"],"");
end    
sgtitle(sgtitle_text,"fontsize",23,"interpreter","latex")

figure(2)
leg = legend(string(2*Ns),'Interpreter','latex','linewidth',2);
title(leg,"$N$")

fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.92;
Lgnd.Position(2) = 0.4;

sgtitle_text = join([curve_name,", $|M|$ = ",len,", $\mu=$ ",num2str(Mu),", 'true' $N$=",2*N2],"");
if rank_def ~=0
sgtitle_text = join([sgtitle_text,", SVD"],"");
end    
sgtitle(sgtitle_text,"fontsize",23,"interpreter","latex")
% [string(2*Ns(1:end-1)); string(2*Ns(2:end))]
