% This function returns opt sig and sig'
% Inputs : curve name, point set = t, normalize length = len1, curve
% params, mu, k = evalue #, V = velocities as cells of 2 rows (1st row dx/da, 2nd row dy/da)
% Outputs : sigk, dsigk
function sigk = sh_bayes_fn(curve_number,curve_params,mu,N1,N,M,len1,tol,Hom,rho_wt,k,mult,pm,method,pn)
    
    [t,x,dx,d2x,nx,points_inside,meshes,~,~,~,dp2norm,len,dom_area,curve_name] = point_setup(curve_number,curve_params,N,M,len1,Hom);
    switch method
        case 0
            [~,evs] = alg1_sh(curve_number,curve_params,mu,N1,N,M,len1,tol,Hom,rho_wt,0);
            sigk = pm*evs(k:(k+mult-1))';
        case 1
            [~,evs] = alg1_sh(curve_number,curve_params,mu/sqrt(dom_area),N1,N,M,len1,tol,Hom,rho_wt,0);
            sigk = pm*evs(k:(k+mult-1))';
            sigk = sigk*len;
        case 2
            [~,evs] = alg1_sh(curve_number,curve_params,mu/len,N1,N,M,len1,tol,Hom,rho_wt,0);
            %sigk = pm*evs(k:(k+mult-1))';
            if pn == "n"
            sigk = pm*max(evs(evs < 0)); 
            else
            sigk = pm*min(evs(evs > 0)); 
            end
            sigk = sigk*sqrt(dom_area);
        case 3
            [~,evs] = alg1_sh(curve_number,curve_params,mu/len,N1,N,M,len1,tol,Hom,rho_wt,0);
            % sigk = pm*evs(k:(k+mult-1))';
            if pn == "n"
            sigk = pm*max(evs(evs < 0)); 
            else
            sigk = pm*min(evs(evs > 0)); 
            end
            sigk = sigk*len;
        case 4
            [~,evs] = alg1_sh(curve_number,curve_params,mu/len,N1,N,M,len1,tol,Hom,rho_wt,0);
            sigk = pm*evs(k:(k+mult-1))';
            sigk = sigk*dom_area;
        case 5
            [~,evs] = alg1_sh(curve_number,curve_params,sqrt(dom_area)*mu/len,N1,N,M,len1,tol,Hom,rho_wt,0);
            sigk = pm*evs(k:(k+mult-1))';
            sigk = sigk*len/sqrt(dom_area);
        case 6
            [~,evs] = alg1_sh(curve_number,curve_params,len*mu/sqrt(dom_area),N1,N,M,len1,tol,Hom,rho_wt,0);
            sigk = pm*evs(k:(k+mult-1))';
            sigk = sigk*sqrt(dom_area)/len;    
    end
    % sigk = pm*evs(k:(k+mult-1))';
    % sigk = sigk*len;
    
    % sigk = pm*evs(k:(k+mult-1))';
    % sigk = sigk*sqrt(dom_area);
end