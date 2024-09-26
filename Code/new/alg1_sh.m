% Algorithm 1 for SH EVP 
function [rank_def,evs,eden,evecs,A,B] = alg1_sh(curve_number,curve_params,Mu,N1,N,M,len1,tol,Hom,rho_wt,eigsk)
[rank_def,~,~,~,~,~] = stek_helm(curve_number,curve_params,Mu,N1,M,len1,tol,1,Hom,rho_wt,eigsk);
if rank_def ~= 0 
    [rank_def,evs,eden,evecs,A,B] = stek_helm(curve_number,curve_params,Mu,N,M,len1,tol,1,Hom,rho_wt,eigsk);
else
    [rank_def,evs,eden,evecs,A,B] = stek_helm(curve_number,curve_params,Mu,N,M,len1,tol,0,Hom,rho_wt,eigsk);
end    

end