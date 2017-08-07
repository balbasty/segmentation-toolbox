function mu = sample_mu(mu,phi)

deg = 2;
d   = size(mu);
Kb  = d(4);

mu = double(mu);
x1 = double(phi(:,:,:,1));
x2 = double(phi(:,:,:,2));
x3 = double(phi(:,:,:,3));

for k1=1:Kb
    bg1(k1) = mean(mean(mu(:,:,1,k1)));
    bg2(k1) = mean(mean(mu(:,:,end,k1)));
end

s    = cell(1,Kb);
msk1 = x1>=1 & x1<=d(1) & x2>=1 & x2<=d(2) & x3>=1 & x3<=d(3);
msk2 = x3<1;
x1 = x1(msk1);
x2 = x2(msk1);
x3 = x3(msk1);

tot = zeros(d(1:3));
for k=1:Kb,
    a    = spm_bsplins(mu(:,:,:,k),x1,x2,x3,[deg deg deg  0 0 0]);
    s{k} = ones(d(1:3))*bg2(k);
    s{k}(msk1) = exp(a);
    s{k}(msk2) = bg1(k);
    tot  = tot + s{k};
end
msk      = ~isfinite(tot);
tot(msk) = 1;
for k=1:Kb,
    s{k}(msk) = bg2(k);
    s{k}      = s{k}./tot;
end

mu = zeros(d,'single');
for k=1:Kb
   mu(:,:,:,k) = s{k};
end