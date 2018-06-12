function ElnPzN = compute_ElnPzN(mrf,lkp,resp)
Kb = max(lkp.part);
R  = zeros([mrf.dm Kb],'single');
for z=1:mrf.dm(3)
    for k=1:Kb
        slice = 0;
        for k1=find(lkp.part==k)
            slice = slice + resp(k1).dat(:,:,z);
        end
            
        R(:,:,z,k) = slice;     
    end   
    clear slice
end

if mrf.ml
    ElnPzN = spm_vbmrf_lowerbound(R,single(log(mrf.Upsilon)),mrf.w);
else
    ElnPzN = spm_vbmrf_lowerbound(R,single(mrf.ElnUpsilon),mrf.w);
end
%==========================================================================