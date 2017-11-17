function mom = mom_John2Bishop(mom)
K = numel(mom.s0);
for k=1:K
    s1(:,k)   = mom.s1(:,k)/mom.s0(k);
    S2(:,:,k) = mom.S2(:,:,k)/mom.s0(k) - (mom.s1(:,k)/mom.s0(k))*(mom.s1(:,k)/mom.s0(k))';
end
mom.s1 = s1;
mom.S2 = S2;
%==========================================================================