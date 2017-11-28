function M = threshold_eig(M)
[V,D1] = eig(M);
tol    = max(diag(D1))*eps('single');
D1     = diag(max(diag(D1),tol));
M      = real(V*D1*V');