syms x a lam K

syms b1 b2 b3 b4
syms c1 c2 c3 c4
syms w1 w2 w3 w4

%%
F = a/x - K*( (b1/(w1*x +c1)) + (b2/(w2*x +c2)) + (b3/(w3*x +c3)) + (b4/(w4*x +c4)) ) - lam;

solve(F,x)

%%
F = b1*lam + sqrt(lam^2*b1^2 + lam*b3)