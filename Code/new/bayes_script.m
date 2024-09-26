run inputs.m

meths = 2:3; exps = 20; 
aa = -1; bb = 1; 
ev_ks = 1;%1:5; 

mus = [1 5 7];
PN = ["n", "p"];
mm = ["max","min"]; 
PMs = [-1 1]; Eigt = ["neg","pos"]; 
Mag = ["largest","smallest"]; % ""

for Mu = 1:length(mus)
mu = mus(Mu);

for pns = 1:2
pn = PN(pns); mag = Mag(pns); eigt=Eigt(pns);
for pms = 1:2
pm = PMs(pms);
maxmin = mm(pms);

for evk = ev_ks
ev_k = evk;    
nvar = 1+length(curve_params); % mu + cp
varnames = cell(1,nvar-1); vars = cell(1,nvar-1);
nvar = nvar-1; % only cp
mu_x_fun = zeros(length(meths),nvar+1);
trial_x_fun = zeros(length(exps),nvar+1);
cp_type = 'real';

var0 = optimizableVariable('mu',[1,20],'Type','real');
for i = 1:numel(varnames)
    varnames{i} = join(["cp",i],"");
    vars{i} = optimizableVariable(varnames{i}, [aa bb], 'Type', cp_type);
end

% mu = 5;
for ms = meths
% fun = @(x)sh_bayes_fn(curve_number,[x.cp1 x.cp2 x.cp3 x.cp4 x.cp5 x.cp6 x.cp7 x.cp8 x.cp9 x.cp10 x.cp11 x.cp12 x.cp13 x.cp14 x.cp15 x.cp16],x.mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms);
% fun = @(x)sh_bayes_fn(curve_number,[x.cp1 x.cp2 x.cp3 x.cp4 x.cp5 x.cp6 x.cp7 x.cp8 x.cp9 x.cp10],x.mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms);
fun = @(x)sh_bayes_fn(curve_number,[x.cp1 x.cp2 x.cp3 x.cp4 x.cp5 x.cp6],mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms,pn);
% fun = @(x)sh_bayes_fn(curve_number,getVariableValues(x,varnames),x.mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms);
% fun = @(x)sh_bayes_fn(curve_number,[x.cp1],x.mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms);
% fun = @(x)sh_bayes_fn(curve_number,[x.cp1],mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms);
for e = 1:exps
%results = bayesopt(fun,[var0,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10]); 
% results = bayesopt(fun,[var0,[vars{:}]]); 
results = bayesopt(fun,[vars{:}]); 
%mu_x_fun(ms+1,:) = [table2array(results.XAtMinObjective) results.MinObjective]; 
trial_x_fun(e,:) = [table2array(results.XAtMinObjective) results.MinObjective]; 
close all
disp(join(["mu = ",mu,", pn = ",pn,", pm = ", pm,", method  = ", ms, ", exp = ",e],""))
end
if exist(join(["opt/files/",curve_name],""),'dir') == 0
    mkdir(join(["opt/files/",curve_name],""))
end    
% if pm == 1
% filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nvar,".csv"],"");
    % filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nvar,"_min_largest_neg.csv"],"");
    filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nvar,"_",maxmin,"_",mag,"_",eigt,".csv"],"");
% else
%     % filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nvar,"_max_largest_neg.csv"],"");
%     filename = join(["opt/files/",curve_name,"/mu_x_fnval_ev",ev_k,"_method",ms,"_trials",exps,"_",curve_name,"_nparams",nvar,"_",maxmin,"_",mag,"_",PM,".csv"],"");
% end    

if nvar == length(curve_params)
    filename = strrep(filename,"_nparams",join(["_mu",mu,"_nparams"],""));
end    

writematrix(trial_x_fun,filename)
end
% writematrix(mu_x_fun,join(["opt/files/mu_x_fnval_ev",ev_k,".csv"],""))
end
end % pm end
end % mu end
end % pn end