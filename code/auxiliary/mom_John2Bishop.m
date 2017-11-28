function mom = mom_John2Bishop(mom)
K = numel(mom(1).s0);
for i=1:numel(mom) 
    s1 = zeros(size(mom(i).s1));
    S2 = zeros(size(mom(i).S2));
    for k=1:K
        s1(:,k)   = mom(i).s1(:,k)/mom(i).s0(k);
        S2(:,:,k) = mom(i).S2(:,:,k)/mom(i).s0(k) - (mom(i).s1(:,k)/mom(i).s0(k))*(mom(i).s1(:,k)/mom(i).s0(k))';
    end
    mom(i).s1 = s1;
    mom(i).S2 = S2;
end
%==========================================================================