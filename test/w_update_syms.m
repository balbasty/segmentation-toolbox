% TPM updates without scaling by mixing proportions
syms  a1  a2  a3  a4
syms  b11 b12 b13 b14
syms  b21 b22 b23 b24
syms  b31 b32 b33 b34

a  = [a1  a2  a3  a4].';
b1 = [b11 b12 b13 b14].';
b2 = [b21 b22 b23 b24].';
b3 = [b31 b32 b33 b34].';

syms  c11 c12 c13 c14
syms  c21 c22 c23 c24
syms  c31 c32 c33 c34
c1 = [c11 c12 c13 c14].';
c2 = [c21 c22 c23 c24].';
c3 = [c31 c32 c33 c34].';

q1 = exp(a+c1)./sum(exp(a+c1));
q2 = exp(a+c2)./sum(exp(a+c2));
q3 = exp(a+c3)./sum(exp(a+c3));
F  = log(sum(b1.*q1)) + log(sum(b2.*q2)) + log(sum(b3.*q3));

dF = simplify(diff(F,c11))