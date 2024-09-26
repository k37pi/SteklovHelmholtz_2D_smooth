run inputs.m

% m0 = nothing(nobody asked what is in your head), 
% m1 = len sig (mu/sq area)
% m2 = sq area sig (mu/len), 
% m3 = len sig (mu/len)
% m4 = area sig (mu/len) [expected for convex domains]
% m5 = len/sqrt(area) sig(mu*sqrt(area)/len)
% m6 = sqrt(area)/len sig(len*mu/sqrt(area))

ev_k = 1; Ms = 2:3; 
exps = 20; 
mu_fixed = 1; % 0 no, 1 yes 
mus = 5;[1 5 7];
Nparams = 10;   2*(1:5);length(curve_params)+1-mu_fixed; 
nrow = exps/5;
scale_text = ["$\sqrt{|\Omega|}\sigma(\Omega,\mu/|\Gamma|)$","${|\Gamma|}\sigma(\Omega,\mu/|\Gamma|)$"];

mm = ["max","min"]; 
PMs = [-1 1]; Eigt = ["neg","pos"];  % "pos" +, "neg" - 
Mag = ["largest","smallest"]; % ""

for NP = 1:length(Nparams) % nparams
nparams = Nparams(NP);
for MM = 1:length(mm) % max / min
maxmin = mm(MM); 

for val = 1:length(Eigt) % pos / neg
eigt = Eigt(val);
mag = Mag(val); % "smallest" for "pos", "largest" for "neg"


multi = 1;
if maxmin == "max"
    multi = -multi;
end    

for Mus = 1:length(mus) % wave number
mu = mus(Mus);
for ms = Ms    % scaling method
%filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nparams,".csv"],"");
if mu_fixed == 1 % mu is fixed
filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_mu",mu,"_nparams",nparams,"_",maxmin,"_",mag,"_",eigt,".csv"],"");
else
 filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nparams,"_",maxmin,"_",mag,"_",eigt,".csv"],"");
end
trial_x_fun = readmatrix(filename);

if mu_fixed == 1
cps = trial_x_fun(:,1:end-1); 
else    
mus = trial_x_fun(:,1);
cps = trial_x_fun(:,2:end-1); 
end
    
funvals = trial_x_fun(:,end);
figure
for c = 1:exps % cth trial
[t,x,dx,d2x,nx,points_inside,meshes,X,R,dp1norm,dp2norm,len,dom_area,curve_name] = point_setup(curve_number,cps(c,:),N,M,len1,Hom);
subplot(nrow,ceil(exps/nrow),c)
if abs(funvals(c)) == abs(max(funvals))
plot(x(1,:),x(2,:),'r'); axis equal; axis tight
elseif abs(funvals(c)) == abs(min(funvals))
plot(x(1,:),x(2,:),'b'); axis equal; axis tight 
else
plot(x(1,:),x(2,:),'k'); axis equal; axis tight     
end
if mu_fixed == 1
title(join(["$\sigma_",ev_k,"$ = ",multi*round(funvals(c),3)],""),FontSize=18)
else    
title(join(["$\sigma_",ev_k,"$ = ",multi*round(funvals(c),3),", $\mu$ = ",round(mus(c),2)],""),FontSize=18)
end
xlabel(join(["$|M|$ = ",round(len,2),", $ \sqrt{|\Omega|} $ = ",round(sqrt(dom_area),2)],""),FontSize=14)
if ms ~= 6
ylabel(join(["$|M|/\sqrt{|\Omega|}$ = ",round(len/sqrt(dom_area),3)],""),FontSize=14)
else
ylabel(join(["$\sqrt{|\Omega|}/|M|$ = ",round(sqrt(dom_area)/len,3)],""),FontSize=14)
end    
end 
if mu_fixed == 1
sgtitle(join(["number of parameters = ",size(cps,2),", scaling = ",scale_text(ms-1),", $\mu$ = ",mu,", ",join([maxmin,mag,eigt]," ")],""),fontsize=20)
else
sgtitle(join(["number of parameters = ",size(cps,2),"+1 (for $\mu$), scaling = ",scale_text(ms-1),", ",join([maxmin,mag,eigt],"_")],""),fontsize=20)
end    

plotname = strrep(filename,".csv",".png");
if exist(join(["opt/files/",curve_name,join(["/nparams",nparams],"")],""),'dir') == 0
    mkdir(join(["opt/files/",curve_name,join(["/nparams",nparams],"")],""))
end    
plotname = strrep(plotname,join(["opt/files/",curve_name],""),join(["opt/files/",curve_name,join(["/nparams",nparams],"")],""));

% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
exportgraphics(gcf,plotname,'Resolution',300) 
saveas(gcf,plotname);
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
hgexport(ax, strrep(plotname,".png",".eps"),...
    hgexport('factorystyle'), 'Format', 'eps');
end
end
end
end
end