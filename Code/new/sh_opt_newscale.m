run inputs.m
ev_k = 1; mult = 1; len1 = 0; 
pm = -1;
MS = 3;
pn = "n";
NP = 2*(1:8);

% minr = 1; maxr = 3; 
for Ms = 1:length(MS)
ms = MS(Ms);
for ns = 1:length(NP)
nparam = NP(ns);length(curve_params);

exps = 1; %exp_deets = zeros(exps,4); exp_deets(:,1) = 1:exps;
% Minr = [0.1,0.1]; % ellipise a min = 3, ellipise b min = 2
% Maxr = [10,10]; % ellipise a max = 8, ellipise b max = 6

Minr = zeros(1,nparam)-.1; % min of each parameter is -0.2;[0.1 0 0 0 0.1 0 0.1 0 0 -0.1]-0.05;zeros(1,nparam)-0.2; % min of each parameter is -0.2
Maxr = zeros(1,nparam)+.1; % min of each parameter is -0.2[0.1 0 0 0 0.1 0 0.1 0 0 -0.1]+0.05;zeros(1,nparam)+0.2; % max of each parameter is 0.2

%Minr = minr; Maxr = maxr;

nonlcon = []; alg = "sqp";

%so_fun = @(cp0) sigopt_new(curve_number,cp0,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm);
its = cell(exps,1); fvals = cell(exps,1); dfvals = cell(exps,1);

%parfor ex = 1:exps
% for ex = 1:exps
% cp = Minr+(Maxr-Minr).*rand(1,nparam);
% cp_fix = cp;
so_fun = @(cp0) sh_bayes_fn(curve_number,cp0,mu,N1,N,M,len1,tol,Hom,rho_wt,ev_k,mult,pm,ms,pn);
% sigopt_new(curve_number,cp0,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm);
psoptions = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','HybridFcn',@fmincon,'display','iter');
[config,fopts,exitflag,op] = particleswarm(so_fun,nparam,Minr,Maxr,psoptions);
filename = join(["opt/files/",curve_name,"/pso_config_method",ms,"_min_neg_",curve_name,"_nparams",nparam,".csv"],"");
filename = strrep(filename,"_nparams",join(["_mu",mu,"_nparams"],""));
writematrix(config,filename);
end
end
% %so_fun = @(cp0) sigopt_new_fixomega(curve_number,cp_fix,cp0,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm);
% 
% disp(cp);
% %[param_opt,fval,maxfval,exitflag,output] = fmincon(so_fun,cp,[],[],[],[],minr,maxr,[],options);
% % exp_deets2(ex) = cp;
% % exp_deets3(ex) = param_opt;
% % exp_deets4(ex) = fval;
% 
% %[xsol,fval,history_fmincon,searchdir] = runfmincon(so_fun,cp,Minr,Maxr,nonlcon,alg);
% [xsol,fval,history_fmincon] = runfmincon(so_fun,cp,Minr,Maxr,nonlcon,alg);
% its{ex} = history_fmincon.x; fvals{ex} = history_fmincon.fval; dfvals{ex} = history_fmincon.dfval;
% disp(ex/exps*100);
% end
% 
% M1 = max(cellfun(@numel, its))/length(cp);
% iterations = cell2mat(its);
% iterations = [repelem(1:exps,cellfun(@numel, its)/length(cp))',reshape(iterations,[],nparam)];
% %iterations = reshape(cell2mat(cellfun(@(x) [x; NaN(M1-numel(x),2)],its,'uni',0)),[],exps);
% 
% M2 = max(cellfun(@numel, fvals));%/length(cp);
% % fevaluations = reshape(cell2mat(cellfun(@(x) [x; NaN(M2-numel(x),1)],fvals,'uni',0)),[],exps);
% fevaluations = cell2mat(fvals);
% fevaluations = [repelem(1:exps,cellfun(@numel, fvals))',fevaluations];
% 
% M3 = max(cellfun(@numel, dfvals));%/length(cp);
% % dfevaluations = reshape(cell2mat(cellfun(@(x) [x; NaN(M3-numel(x),1)],dfvals,'uni',0)),[],exps);
% dfevaluations = cell2mat(dfvals);
% dfevaluations = [repelem(1:exps,cellfun(@numel, dfvals))',dfevaluations];

% if curve_name == "disk_def"
% fname = "Deformed ellipse";
% param = "$\kappa$";
% else
% fname = curve_name;    
% param = "$R$";
% end    
% 
% %if 1==2
% figure(); 
% subplot(1,3,1)
% plot(0:(M1-1),iterations,'--o','LineWidth',2)
% xlabel("$n$",'Interpreter','latex')
% ylabel(param,'Interpreter','latex')
% set(gca,"FontSize",20)
% 
% subplot(1,3,2)
% plot(0:(M2-1),fevaluations,'--o','LineWidth',2)
% xlabel("$n$",'Interpreter','latex')
% ylabel("$\sigma_k$",'Interpreter','latex')
% set(gca,"FontSize",20)
% 
% subplot(1,3,3)
% plot(0:(M3-1),dfevaluations,'--o','LineWidth',2)
% xlabel("$n$",'Interpreter','latex')
% ylabel("$\sigma'_k$",'Interpreter','latex')
% set(gca,"FontSize",20)
% 
% sgtitle(join(["Iterations of ",exps," experiments for ",fname,"$, \mu=$ ",mu,"$, k = $",k,", Norm $|M|=$ ",len1],""),'interpreter','latex',"FontSize",20)
% filename = ['/home/kshitij/Desktop/SFU/KP Thesis draft sfu format/figs/opt/'];
% filename = join([filename,"expts_",curve_name,"_opt_mu",mu,"_k",k,"_len_",len1,"_N",N],'');
% filename = strrep(filename,'.',"dot");
% %filename = join([filename,"_nosecondterm.png"],"");
% filename = join([filename,".png"],"");
% if isfile(filename)
% filename = strrep(filename,'.png',join(["_",randi(1),".png"],""));
% end    
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 25]); 
% 
% for j = 1:5
% k=j;
% cs = 0:1e-1:0.9;-0.5:1e-1:0.5;%0:1e-1:0.9;
% s = zeros(1,length(cs)); ds = s; lens = zeros(1,length(cs));
% s1 = zeros(1,length(cs)); ds1 = s; lens1 = s;
% 
% cs0 = 3;
% for i = 1:length(cs)
%     [~,~,~,~,~,~,~,~,~,~,~,lens(i)] = point_setup(curve_number,cs(i),N,M,len1,Hom);
%     [s(i),ds(i)] = sigopt_new(curve_number,cs(i),mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm);
%     %[s1(i),ds1(i)] = sigopt_new_fixomega(curve_number,cs0,cs(i),mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm);
% end    
% %figure
% %subplot(1,5,j)
% plot(cs,s,'--*'); hold on
% % plot(cs,ds,'--*');
% % plot(cs,lens,'--*')
% % plot(cs(1:end-1),diff(s)./(diff(cs)),'--*');
% % 
% % subplot(1,2,2)
% % plot(cs,s1,'--*'); hold on
% % plot(cs,ds1,'--*');
% % plot(cs(1:end-1),diff(s1)./(diff(cs)),'--*');
% % 
% % subplot(1,3,3)
% % plot(cs,s,'--*'); hold on
% % plot(cs,s1,'--*');
% 
% disp(j)
% end
% plot(cs,ds,'--*');
% plot(cs,lens,'--*')
% leg_text = string([1:5,"$|M|'$","$|M|$"]);
% legend(leg_text,'interpreter','latex')
% cn = strrep(curve_name,"_","\_");
% title_text = join([cn,", $\mu\in$ [",min(mu./lens),",",max(mu./lens),"]"],"");
% title(title_text,'Interpreter','latex','FontSize',25)
% xlabel("domain parameter")
% 
% alpha0 = 1; rho = 0.9; c = 1e-4;
% xk = cs(5);
% fk = s(5); dfk = ds(5); pk = -dfk;
% alpha = alpha0;
% xplus = xk+alpha*pk; 
% fkplus = sigopt_new(curve_number,xplus,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm); 
% 
% while fkplus > fk + c*alpha*dfk'*pk
%     alpha = rho*alpha;
%     xplus = xk+alpha*pk; 
%     fkplus = sigopt_new(curve_number,xplus,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm); 
% end    
% 
% while abs(fkplus - fk) > 1e-3
% dfk = (fkplus - fk)/(xplus-xk);
% fk = fkplus; 
% while fkplus > fk + c*alpha*dfk'*pk
%     alpha = rho*alpha;
%     xplus = xk+alpha*pk; 
%     fkplus = sigopt_new(curve_number,xplus,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm); 
% end    
% 
% end    