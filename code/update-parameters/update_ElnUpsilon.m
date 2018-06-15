function [mrf,ll,L,resp] = update_ElnUpsilon(mrf,lkp,resp,llr,llrb,buf,mg,gmm,wp,wp_l,fig,L,print_ll)
Kb = max(lkp.part);

Upsalpha = zeros(Kb);
for k=1:Kb
    for l=1:Kb
        Upsalpha(k,l) = spm_vbmrf_update_Upsilon(resp.dat,mrf.w,k,l);
    end
end
clear R

if mrf.ml
    % Update Upsilon
    %----------------------------------------------------------------------
    Upsalpha    = max(Upsalpha,eps);
    mrf.Upsilon = bsxfun(@rdivide,Upsalpha,sum(Upsalpha,2));   
else
    % Update ElnUpsilon
    %----------------------------------------------------------------------
    
    Upsalpha = Upsalpha + mrf.Upsalpha0;
    % EUpsilon = bsxfun(@rdivide,Upsalpha,sum(Upsalpha,2));                          
    
    mrf.ElnUpsilon = zeros(Kb);
    for k=1:Kb
        mrf.ElnUpsilon(k,:) = psi(Upsalpha(k,:)) - psi(sum(Upsalpha(k,:)));
    end    
end

% Update lower bound
%--------------------------------------------------------------------------
ll = llr + llrb;

[mom,dll,mrf,~,resp] = compute_moments(buf,lkp,mg,gmm,wp,wp_l,resp,mrf);        
ll                   = ll + dll;         

dll = spm_VBGaussiansFromSuffStats(mom,gmm);
ll  = ll + sum(sum(dll));  

my_fprintf('Upsilon:%g\t%g\t%g\n',ll,llr,llrb,print_ll);
L{1}(end + 1) = ll;

debug_view('convergence',fig{4},lkp,buf,L);
%==========================================================================