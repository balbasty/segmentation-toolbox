%% My attempt
syms  mu1  mu2  mu3  mu4  % Tissue priors at some voxel

syms  r11 r12 r13 r14 % responsibility for 1st subject
syms  r21 r22 r23 r24 % responsibility for 2nd subject
syms  r31 r32 r33 r34 % responsibility for 3rd subject

syms  w11 w12 w13 w14 % "Mixing weights" for the different subjects
syms  w21 w22 w23 w24
syms  w31 w32 w33 w34

mu  = [mu1  mu2  mu3  mu4].';

r1 = [r11 r12 r13 r14].';
r2 = [r21 r22 r23 r24].';
r3 = [r31 r32 r33 r34].';

w1 = [w11 w12 w13 w14].';
w2 = [w21 w22 w23 w24].';
w3 = [w31 w32 w33 w34].';

q1 = exp(mu+w1)./sum(exp(mu+w1));
q2 = exp(mu+w2)./sum(exp(mu+w2));
q3 = exp(mu+w3)./sum(exp(mu+w3));

F = log(sum(r1.*q1)) + log(sum(r2.*q2)) + log(sum(r3.*q3));

% o1=w1(1)/sum(mu.*w1);
% o2=w2(1)/sum(mu.*w2);
% o3=w3(1)/sum(mu.*w3);

p1=1/sum(exp(mu + w1));
p2=1/sum(exp(mu + w2));
p3=1/sum(exp(mu + w3));

a1=r1.*exp(mu+w1)/sum(r1.*exp(mu+w1));
a2=r2.*exp(mu+w2)/sum(r2.*exp(mu+w2));
a3=r3.*exp(mu+w3)/sum(r3.*exp(mu+w3));

%% Check that the following comes out as zero...
simplify(a1(1) - exp(mu1 + w11)*p1  +  a2(1) - exp(mu1 + w21)*p2  +  a3(1) - exp(mu1 + w31)*p3  -  diff(F,mu1))

%% Updates via...
solve(sym('a1(1) - exp(mu1 + w11)*p1  +  a2(1) - exp(mu1 + w21)*p2  +  a3(1) - exp(mu1 + w31)*p3'),mu1)