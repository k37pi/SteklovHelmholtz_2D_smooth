%% Mus and inputs --- 
run inputs.m
Mus = [0.1:1e-1:20];
Sigs = zeros(2*N,length(Mus)+1);

%% spectra -----
Sigs(:,1) = stek_lap(curve_number,curve_params,N,M,len1,Hom,rho_wt,intext);

for i = 1:length(Mus)
[~,sigs] = alg1_sh(curve_number,curve_params,Mus(i),N1,N,M,len1,tol,Hom,rho_wt,0);
Sigs(1:length(sigs),i+1) = sigs; 
disp(i/length(Mus))
end

%% plot
plot(zeros(N,1),Sigs(1:N,1),'--'); hold on
for i = 1:length(Mus)
plot(Mus(i)+zeros(N,1),Sigs(1:N,i+1),'--o')
disp(i/length(Mus))
end

