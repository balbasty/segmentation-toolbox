function [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,iter,tol1,nm,nitgmm,fig,L,print_ll)
tiny = eps*eps;
K    = numel(mg);
Kb   = numel(wp);

for subit=1:nitgmm
    oll = ll;
    ll  = llrb + llr;
    
    % Compute responsibilities and moments
    [mom,dll,mgm] = compute_moments(buf,K,mg,gmm,wp,lkp);        
    ll            = ll + dll; 
    
    % Compute missing data and VB components of ll
    dll = spm_VBGaussiansFromSuffStats(mom,gmm);
    ll  = ll + sum(sum(dll));
    
    L{1}(end + 1) = ll;
    
    % Add up 0:th moment
    s0 = 0;
    for i=2:numel(mom), s0 = s0 + mom(i).s0; end
    
    % Update tissue weights
    for k1=1:Kb
        wp(k1) = (sum(s0(lkp==k1)) + wp_reg*1)/(mgm(k1) + wp_reg*Kb); % bias the solution towards 1
    end
    wp = wp/sum(wp);
   
    % Update mixing proportions
    for k=1:K
        tmp   = s0(lkp==lkp(k));
        mg(k) = (s0(k) + tiny)/sum(tmp + tiny);  % US eq. 27 (partly)
    end
        
    % Update means and variances
    [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);

    my_fprintf('MOG:\t%g\t%g\t%g\n',ll,llr,llrb,print_ll);

    if subit>1 || iter>1
        debug_view('convergence',fig{4},lkp,buf,L);
    end
    if subit>1 && ll-oll<tol1*nm
        % Improvement is small, so go to next step
        break;
    end
end
%=======================================================================