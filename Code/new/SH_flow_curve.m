%% Curve params and inputs --- 
run inputs.m
cps = 3;0:.5*1e-1:.99;1:10;0:1e-1:1.5;
Mus = 1:20; 
[C,WN] = meshgrid(cps,Mus);
Sigs_dom = zeros(2*N,length(cps));
Sigs_len = zeros(2*N,length(cps));
ks = 20;
f_evlyn = @(cp,Mu) ev_len(curve_number,cp,Mu,N1,N,M,len1,tol,Hom,rho_wt,pm,ks,0);    
tic
[Ss,Ss_D,Ss_LD,Ss_L,Ls,As] = arrayfun(@(x1,x2) f_evlyn(x1,x2), C, WN,'UniformOutput', false);
toc

Lens = [Ls{:}]';
Lens = reshape(Lens,length(Mus),length(cps));

Areas = [As{:}]';
Areas = reshape(Areas,length(Mus),length(cps));

Sig = [Ss{:}]';
Sig = reshape(Sig,length(Mus),length(cps),ks);

SigD = [Ss_D{:}]';
SigD = reshape(SigD,length(Mus),length(cps),ks);

SigLD = [Ss_LD{:}]';
SigLD = reshape(SigLD,length(Mus),length(cps),ks);

SigL = [Ss_L{:}]';
SigL = reshape(SigL,length(Mus),length(cps),ks);

%% -----
cn = strrep(curve_name,"_","\_");"kite";
nn = 1;
for nn = 1:ks
subplot(1,4,1)
surf(C,WN,Sig(:,:,nn)); 
% shading interp
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = join(["$\sigma$",nn,"$(\Omega,\mu)$"],"");
title(title_text);

subplot(1,4,2)
surf(C,WN,SigD(:,:,nn)); 
%shading interp
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = join(["$|M|\sigma",nn,"(\Omega,\mu/\sqrt{|\Omega|})$"],"");
title(title_text);

subplot(1,4,3)
surf(C,WN,SigLD(:,:,nn)); 
%shading interp
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = join(["$\sqrt{|\Omega|}\sigma",nn,"(\Omega,\mu/|M|)$"],"");
title(title_text);

subplot(1,4,4)
surf(C,WN,SigL(:,:,nn)); 
%shading interp
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = join(["$|M|\sigma",nn,"(\Omega,\mu/|M|)$"],"");
title(title_text);

sgtitle(cn)

filename = join(["/home/kshitij/Desktop/SFU/Steklov/opt/",curve_name," grid search/ev",nn,curve_name],"");
ax = gcf;
saveas(gcf,filename);
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
hgexport(ax, filename,...
    hgexport('factorystyle'), 'Format', 'eps');
end

%% -----------------
for nn = 1:6
subplot(2,3,nn)
surf(C,WN,SigL(:,:,nn)); 
%shading interp
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = join(["$|M|\sigma",nn,"(\Omega,\mu/|M|)$"],"");
title(title_text);

sgtitle(cn)
end
filename = join(["/home/kshitij/Desktop/SFU/Steklov/opt/",cn," grid search/evL",nn,curve_name],"");
ax = gcf;
saveas(gcf,filename);
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
hgexport(ax, filename,...
    hgexport('factorystyle'), 'Format', 'eps');


figure
subplot(2,2,1)
surf(C,WN,Lens);
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = "$|M|$";
title(title_text);

subplot(2,2,2)
surf(C,WN,sqrt(Areas));
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = "$\sqrt{|\Omega|}$";
title(title_text);

subplot(2,2,3)
surf(C,WN,WN./Lens);
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = "$\mu/|M|$";
title(title_text);

subplot(2,2,4)
surf(C,WN,WN./sqrt(Areas));
colorbar(); colormap("jet")
view(0,90)
xlabel("$c$"); ylabel("$\mu$"); 
title_text = "$\mu/\sqrt{|\Omega|}$";
title(title_text);
filename = join(["/home/kshitij/Desktop/SFU/Steklov/opt/",cn," grid search/params",curve_name],"");
ax = gcf;
saveas(gcf,filename);
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
hgexport(ax, filename,...
    hgexport('factorystyle'), 'Format', 'eps');