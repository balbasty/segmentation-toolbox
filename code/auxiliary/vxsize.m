function vx = vxsize(M)
M  = M(1:3,1:3);
vx = sqrt(sum(M.^2));