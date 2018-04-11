function mg = update_mg(lkp,s0)
tiny = eps*eps;
K    = numel(lkp.part);
mg   = zeros(1,K);
for k=1:K
    tmp   = s0(lkp.part==lkp.part(k));
    mg(k) = (s0(k) + tiny)/sum(tmp + tiny); % US eq. 27 (partly)
end

