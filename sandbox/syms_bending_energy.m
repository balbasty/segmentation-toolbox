% MATLAB code (requiring Symbolic Toolbox) for computing the
% regulariser for bending energy. Note that an additional multiplication
% by the square of the voxel size is needed.

clear; clc;

syms v0 v1 v2

syms s1 s2 s3
syms l1 l2 l3
zz = sym(zeros(3,3));
K1 = cat(3,zz,[0 -1/(v0*v0) 0; 0 2/(v0*v0) 0; 0 -1/(v0*v0) 0],zz);
K2 = cat(3,zz,[0 0 0; -1/(v1*v1) 2/(v1*v1) -1/(v1*v1); 0 0 0],zz);
K3 = sym(zeros(3,3,3));
K3(2,2,1) = -1/(v2*v2);
K3(2,2,2) = 2/(v2*v2);
K3(2,2,3) = -1/(v2*v2);

K  = K1+K2+K3;
K1 = K*l1; K1(2,2,2) = K1(2,2,2)+l2;
K2 = K*l1; K2(2,2,2) = K2(2,2,2)+l3;

% L  = convn(K,K)
L  = sym(zeros(5,5,5));
for i=1:3,
    for j=1:3,
        for k=1:3,
            L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) + K1(i,j,k)*K2;
        end;
    end;
end;
disp(simplify(L(3:end,3:end,3:end)))