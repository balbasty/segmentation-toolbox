function lnPzN = compute_lnPzN(mrf,resp)
% Compute log-MRF part of responsibilities
if mrf.ml
    lnPzN = spm_vbmrf(resp.dat,single(log(mrf.Upsilon)),mrf.w);
else
    lnPzN = spm_vbmrf(resp.dat,single(mrf.ElnUpsilon),mrf.w);
end
%==========================================================================