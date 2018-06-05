function mrf = init_mrf(obj,d,lkp,vx)
K      = numel(lkp.part);
Kb     = max(lkp.part);
mrf    = obj.segment.mrf;
mrf.dm = d;    

if ~isfield(mrf,'ElnPzN')
    mrf.ElnPzN = 0;
end

if isempty(mrf.lnUpsilon)
    val_diag = mrf.val_diag;        
    Upsilon  = ones(Kb,'single');
    Upsilon  = Upsilon + (val_diag - 1)*eye(Kb);
    Upsilon  = bsxfun(@rdivide,Upsilon,sum(Upsilon,2));

    mrf.lnUpsilon = log(Upsilon);
end

mrf.vx = single(1./(vx.^2));
%==========================================================================