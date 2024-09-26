function [evs,evs_d,evs_ld,evs_l,len,dom_area] = ev_len(curve_number,curve_params,Mu,N1,N,M,len1,tol,Hom,rho_wt,pm,k,eigsk)
[~,~,~,~,~,~,~,~,~,~,~,len,dom_area,~] = point_setup(curve_number,curve_params,N,M,len1,Hom);

[~,evs] = alg1_sh(curve_number,curve_params,Mu,N1,N,M,len1,tol,Hom,rho_wt,eigsk);
evs = pm*evs(1:k);

[~,evs_d] = alg1_sh(curve_number,curve_params,Mu/sqrt(dom_area),N1,N,M,len1,tol,Hom,rho_wt,eigsk);
evs_d = pm*evs_d(1:k)*len;

[~,evs_ld] = alg1_sh(curve_number,curve_params,Mu/len,N1,N,M,len1,tol,Hom,rho_wt,eigsk);
evs_ld = pm*evs_ld(1:k)*sqrt(dom_area);

[~,evs_l] = alg1_sh(curve_number,curve_params,Mu/len,N1,N,M,len1,tol,Hom,rho_wt,eigsk);
evs_l = pm*evs_l(1:k)*len;
end

