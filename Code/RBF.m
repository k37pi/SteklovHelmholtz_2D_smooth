function rbfx = RBF(points,rbf,rbf_coeffs,rbf_points)
    
    rbfx = zeros(length(points),1);
    for i = 1:length(points)
    point = points(i); 
    for j = 1:length(rbf_points)
        distr = norm(point-rbf_points(j)); 
        rbfx(i) = rbfx(i)+rbf_coeffs(j)*rbf(distr);
    end    
    end

end