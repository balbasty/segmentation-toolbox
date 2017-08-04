%% My attempt
syms  lnmu1  lnmu2  lnmu3  lnmu4  % logs of tissue priors at some voxel

syms  r11 r12 r13 r14 % responsibility for 1st subject
syms  r21 r22 r23 r24 % responsibility for 2nd subject
syms  r31 r32 r33 r34 % responsibility for 3rd subject

syms  lnw11 lnw12 lnw13 lnw14 % logs of "Mixing weights" for the different subjects
syms  lnw21 lnw22 lnw23 lnw24
syms  lnw31 lnw32 lnw33 lnw34

lnmu  = [lnmu1  lnmu2  lnmu3  lnmu4].';

r1 = [r11 r12 r13 r14].';
r2 = [r21 r22 r23 r24].';
r3 = [r31 r32 r33 r34].';

lnw1 = [lnw11 lnw12 lnw13 lnw14].';
lnw2 = [lnw21 lnw22 lnw23 lnw24].';
lnw3 = [lnw31 lnw32 lnw33 lnw34].';

q1 = exp(lnmu + lnw1)./sum(exp(lnmu + lnw1));
q2 = exp(lnmu + lnw2)./sum(exp(lnmu + lnw2));
q3 = exp(lnmu + lnw3)./sum(exp(lnmu + lnw3));

F = sum(r1.*log(q1)) + sum(r2.*log(q2)) + sum(r3.*log(q3));

dF = diff(F,lnmu1);

sdF = simplify(dF);

% solve(dF,lnmu1)

%% r_{mk}
a = (r11*(exp(lnmu2 + lnw12) + exp(lnmu3 + lnw13) + exp(lnmu4 + lnw14)))/(exp(lnmu1 + lnw11) + exp(lnmu2 + lnw12) + exp(lnmu3 + lnw13) + exp(lnmu4 + lnw14)) ...
    + (r21*(exp(lnmu2 + lnw22) + exp(lnmu3 + lnw23) + exp(lnmu4 + lnw24)))/(exp(lnmu1 + lnw21) + exp(lnmu2 + lnw22) + exp(lnmu3 + lnw23) + exp(lnmu4 + lnw24)) ...
    + (r31*(exp(lnmu2 + lnw32) + exp(lnmu3 + lnw33) + exp(lnmu4 + lnw34)))/(exp(lnmu1 + lnw31) + exp(lnmu2 + lnw32) + exp(lnmu3 + lnw33) + exp(lnmu4 + lnw34)) ...
    - (r12*exp(lnmu1 + lnw11))/(exp(lnmu1 + lnw11) + exp(lnmu2 + lnw12) + exp(lnmu3 + lnw13) + exp(lnmu4 + lnw14)) ...
    - (r13*exp(lnmu1 + lnw11))/(exp(lnmu1 + lnw11) + exp(lnmu2 + lnw12) + exp(lnmu3 + lnw13) + exp(lnmu4 + lnw14)) ...
    - (r14*exp(lnmu1 + lnw11))/(exp(lnmu1 + lnw11) + exp(lnmu2 + lnw12) + exp(lnmu3 + lnw13) + exp(lnmu4 + lnw14)) ...
    - (r22*exp(lnmu1 + lnw21))/(exp(lnmu1 + lnw21) + exp(lnmu2 + lnw22) + exp(lnmu3 + lnw23) + exp(lnmu4 + lnw24)) ...
    - (r23*exp(lnmu1 + lnw21))/(exp(lnmu1 + lnw21) + exp(lnmu2 + lnw22) + exp(lnmu3 + lnw23) + exp(lnmu4 + lnw24)) ...
    - (r24*exp(lnmu1 + lnw21))/(exp(lnmu1 + lnw21) + exp(lnmu2 + lnw22) + exp(lnmu3 + lnw23) + exp(lnmu4 + lnw24)) ...
    - (r32*exp(lnmu1 + lnw31))/(exp(lnmu1 + lnw31) + exp(lnmu2 + lnw32) + exp(lnmu3 + lnw33) + exp(lnmu4 + lnw34)) ...
    - (r33*exp(lnmu1 + lnw31))/(exp(lnmu1 + lnw31) + exp(lnmu2 + lnw32) + exp(lnmu3 + lnw33) + exp(lnmu4 + lnw34)) ...
    - (r34*exp(lnmu1 + lnw31))/(exp(lnmu1 + lnw31) + exp(lnmu2 + lnw32) + exp(lnmu3 + lnw33) + exp(lnmu4 + lnw34));

simplify(a - dF)

%%
b1 = sum(exp(lnmu + lnw1));
b2 = sum(exp(lnmu + lnw2));
b3 = sum(exp(lnmu + lnw3));

na = ( r11*(exp(lnmu2 + lnw12) + exp(lnmu3 + lnw13) + exp(lnmu4 + lnw14)) - (r12 + r13 + r14)*exp(lnmu1 + lnw11) ) / b1 ...
   + ( r21*(exp(lnmu2 + lnw22) + exp(lnmu3 + lnw23) + exp(lnmu4 + lnw24)) - (r22 + r23 + r24)*exp(lnmu1 + lnw21) ) / b2 ...
   + ( r31*(exp(lnmu2 + lnw32) + exp(lnmu3 + lnw33) + exp(lnmu4 + lnw34)) - (r32 + r33 + r34)*exp(lnmu1 + lnw31) ) / b3;

simplify(na - dF)

solve(na,lnmu1)