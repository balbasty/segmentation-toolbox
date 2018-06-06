function mrf = init_mrf(obj,d,lkp,vx)
K      = numel(lkp.part);
Kb     = max(lkp.part);
mrf    = obj.segment.mrf;
mrf.dm = d;    

if ~isfield(mrf,'ElnPzN')
    mrf.ElnPzN = 0;
end

if isempty(mrf.ElnUpsilon)
    val_diag = mrf.val_diag;  
    lambda   = mrf.lambda;
            
    Upsalpha0     = ones(Kb);
    Upsalpha0     = Upsalpha0 + (val_diag - 1)*eye(Kb);
    Upsalpha0     = lambda*Upsalpha0;
    mrf.Upsalpha0 = Upsalpha0;
    
    mrf.ElnUpsilon = zeros(Kb);
    for k=1:Kb
        mrf.ElnUpsilon(k,:) = psi(Upsalpha0(k,:)) - psi(sum(Upsalpha0(k,:)));
    end
end

if ~isfield(mrf,'w')
    mrf.w = single(1./(vx.^2));
end
%==========================================================================