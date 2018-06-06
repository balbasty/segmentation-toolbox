function [mrf,ll,L] = update_ElnUpsilon(mrf,lkp,resp,llr,llrb,buf,mg,gmm,wp,wp_l,fig,L,print_ll)
Kb = max(lkp.part);

% Collect 4D responsibilities
%--------------------------------------------------------------------------
R = zeros([mrf.dm Kb],'single');
for z=1:mrf.dm(3)
    for k=1:Kb
        slice = 0;
        for k1=find(lkp.part==k)
            slice = slice + resp.current(k1).dat(:,:,z);
        end

        R(:,:,z,k) = slice;     
    end    
end
clear slice

Upsalpha = zeros(Kb);
for k=1:Kb
    for l=1:Kb
        Upsalpha(k,l) = spm_vbmrf_update_Upsilon(R,mrf.w,k,l);
    end
end
clear R

Upsalpha = Upsalpha + mrf.Upsalpha0;
% EUpsilon = bsxfun(@rdivide,Upsalpha,sum(Upsalpha,2));                          

% Update ElnUpsilon
%--------------------------------------------------------------------------
mrf.ElnUpsilon = zeros(Kb);
for k=1:Kb
    mrf.ElnUpsilon(k,:) = psi(Upsalpha(k,:)) - psi(sum(Upsalpha(k,:)));
end

% Update lower bound
%--------------------------------------------------------------------------
ll = llr + llrb;

[mom,dll,mrf] = compute_moments(buf,lkp,mg,gmm,wp,wp_l,resp.current,resp.current,mrf);        
ll            = ll + dll;         

dll = spm_VBGaussiansFromSuffStats(mom,gmm);
ll  = ll + sum(sum(dll));  

my_fprintf('Upsilon:%g\t%g\t%g\n',ll,llr,llrb,print_ll);
L{1}(end + 1) = ll;

debug_view('convergence',fig{4},lkp,buf,L);
%==========================================================================