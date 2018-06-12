function mrf = init_mrf(obj,d,lkp,vx)
K      = numel(lkp.part);
Kb     = max(lkp.part);
mrf    = obj.segment.mrf;
mrf.dm = d;    

if ~isfield(mrf,'ElnPzN')
    mrf.ElnPzN = 0;
end

if isempty(mrf.ElnUpsilon) || isempty(mrf.Upsilon)
    if isfield(mrf,'KK')
        KK        = mrf.KK;
        KK(KK==0) = 1e-4;
    end
    
    if mrf.ml
        val_diag = mrf.val_diag;
        
        mrf.Upsilon = ones(Kb);
        mrf.Upsilon = mrf.Upsilon + val_diag*eye(Kb);
        if isfield(mrf,'KK')
            mrf.Upsilon = KK.*mrf.Upsilon;
        end
        
        for k=1:Kb
            mrf.Upsilon(k,:) = mrf.Upsilon(k,:)./sum(mrf.Upsilon(k,:));
        end
    else
        alpha         = mrf.alpha;              
        Upsalpha0     = alpha*ones(Kb);
        mrf.Upsalpha0 = Upsalpha0;

        mrf.ElnUpsilon = zeros(Kb);
        for k=1:Kb
            mrf.ElnUpsilon(k,:) = psi(Upsalpha0(k,:)) - psi(sum(Upsalpha0(k,:)));
        end
    end
end

if ~isfield(mrf,'w')
    mrf.w = single(1./(vx.^2));
end
%==========================================================================