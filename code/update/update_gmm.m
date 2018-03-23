function [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,mix_wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll,wp_lab)
tiny = eps*eps;
K    = numel(mg);
Kb   = numel(wp);

for subit=1:nitgmm
    oll = ll;
    ll  = llrb + llr;
    
    % Compute responsibilities and moments
    [mom,dll,mgm] = compute_moments(buf,K,mg,gmm,wp,lkp,wp_lab);        
    ll            = ll + dll;     
    
    % Add up 0:th moment
    s0 = 0;
    for i=2:numel(mom), s0 = s0 + mom(i).s0; end
    
    nvox = 0;
    for z=1:numel(buf)
        nvox = nvox + buf(z).Nm;
    end
    
    if do_wp
        % Update tissue weights
        for k1=1:Kb
            w1 = mix_wp_reg;
            w2 = 1 - w1;
            
            wp(k1) = (w1*sum(s0(lkp.part==k1)) + w2*nvox*wp_reg)/(w1*mgm(k1) + w2*nvox); % bias the solution towards 1
        end
        wp = wp/sum(wp);
    end
    
    % Update mixing proportions
    for k=1:K
        tmp   = s0(lkp.part==lkp.part(k));
        mg(k) = (s0(k) + tiny)/sum(tmp + tiny);  % US eq. 27 (partly)
    end
        
    % Update means and variances
    [dll,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);
    ll        = ll + sum(sum(dll));  

    my_fprintf('MOG:\t%g\t%g\t%g\n',ll,llr,llrb,print_ll);
    L{1}(end + 1) = ll;

    if subit>1 || iter>1
        debug_view('convergence',fig{4},lkp,buf,L);
    end
    if subit>1 && ll-oll<tol1*nm
        % Improvement is small, so go to next step
        break;
    end
end
%=======================================================================