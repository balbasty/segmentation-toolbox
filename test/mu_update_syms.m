% TPM updates without scaling by mixing proportions
syms  a1  a2  a3  a4
syms  b11 b12 b13 b14
syms  b21 b22 b23 b24
syms  b31 b32 b33 b34

a  = [a1  a2  a3  a4].';
b1 = [b11 b12 b13 b14].';
b2 = [b21 b22 b23 b24].';
b3 = [b31 b32 b33 b34].';

q  = exp(a)./sum(exp(a));
F  = log(sum(b1.*q)) + log(sum(b2.*q)) + log(sum(b3.*q));

p  = 1/sum(exp(a));
r1 = b1.*exp(a)/sum(b1.*exp(a));
r2 = b2.*exp(a)/sum(b2.*exp(a));
r3 = b3.*exp(a)/sum(b3.*exp(a));

% Check that the following comes out as zero...
simplify(r1(1) - exp(a1)*p  +  r2(1) - exp(a1)*p  +  r3(1) - exp(a1)*p  -  diff(F,a1))

% Updates via...
solve(sym('r1(1) - exp(a1)*p  +  r2(1) - exp(a1)*p  +  r3(1) - exp(a1)*p'),a1)


%% TPM updates with scaling by mixing proportions
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

p1=1/sum(exp(a+c1));
p2=1/sum(exp(a+c2));
p3=1/sum(exp(a+c3));
r1=b1.*exp(a+c1)/sum(b1.*exp(a+c1));
r2=b2.*exp(a+c2)/sum(b2.*exp(a+c2));
r3=b3.*exp(a+c3)/sum(b3.*exp(a+c3));

% Check that the following comes out as zero...
simplify(r1(1) - exp(a1 + c11)*p1  +  r2(1) - exp(a1 + c21)*p2  +  r3(1) - exp(a1 + c31)*p3  -  diff(F,a1))

% Updates via...
solve(sym('r1(1) - exp(a1 + c11)*p1  +  r2(1) - exp(a1 + c21)*p2  +  r3(1) - exp(a1 + c31)*p3'),a1)