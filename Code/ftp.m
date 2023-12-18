% This function computes the trig poly a given function f
% Inputs : f, N (# of points = 2*N+1), cos_vecs, sin_vecs 
% cos_vecs, sin_vecs are computed in 'Trig. poly.' of 'main.m' 
% Outputs : tp = trig poly of f evaluated at grid points.

function tp = ftp(f,N,cos_vecs,sin_vecs)
    al = f(:,1:(end-1))*cos_vecs(1:end-1,1:N+1)/N;
    al(:,[1,end]) = al(:,[1,end])/2;
    bl = f(:,1:(end-1))*sin_vecs(1:end-1,2:N)/N;
    tp = al*cos_vecs(:,1:N+1)' + bl*sin_vecs(:,2:N)';
end