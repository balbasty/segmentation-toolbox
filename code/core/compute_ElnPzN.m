function ElnPzN = compute_ElnPzN(mrf,resp)
if mrf.ml
    ElnPzN = spm_vbmrf_lowerbound(resp.dat,single(log(mrf.Upsilon)),mrf.w);
else
    ElnPzN = spm_vbmrf_lowerbound(resp.dat,single(mrf.ElnUpsilon),mrf.w);
end
%==========================================================================