function lnPzN = compute_lnPzN(mrf,lkp,resp)
% Compute log-MRF part of responsibilities
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

lnPzN = spm_vbmrf(R,single(mrf.ElnUpsilon),mrf.w);
%==========================================================================