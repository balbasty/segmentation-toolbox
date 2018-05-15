function obj = modify_bf_dc(obj,iter,verbose)
if nargin<3, verbose = true; end
 
M = numel(obj);    
for m=1:M
    
    if ~obj{m}{1}.segment.do_bf0 || strcmp(obj{m}{1}.modality,'CT')
        continue
    end
    
    S      = numel(obj{m});    
    N      = numel(obj{m}{1}.segment.bf.dc);
    sum_dc = 0;
    for s=1:S    
        % Sum bias field DC components
        sum_dc = sum_dc + obj{m}{s}.segment.bf.dc;
    end
    
    % Average of bias field DC components
    avg_dc = sum_dc./S;
         
    % Mean-correct DC components of bias field
    for s=1:S
        for n=1:N
            obj{m}{s}.segment.Tbias{n}(1,1,1) = obj{m}{s}.segment.Tbias{n}(1,1,1) - avg_dc(n); 
        end
    end
    
    % Convert from log to intensity
    b1 = obj{m}{s}.segment.bf.b1;
    b2 = obj{m}{s}.segment.bf.b2;
    b3 = obj{m}{s}.segment.bf.b3;
 
    scl = exp(b1.*b2.*b3.*avg_dc);
    
    if verbose
        fprintf('%2d | bf_scl = [',iter);
        for n=1:N
            fprintf('%4.2f ',scl(n));
        end
        fprintf(']\n');
    end
    
    % Recompute posteriors with new bias field
    for s=1:S    
        mom = obj{m}{s}.segment.mom;
        gmm = obj{m}{s}.segment.gmm;
        K   = numel(gmm.po.n);
        
        for i=2:numel(mom)
            ind       = mom(i).ind;                        
            scl1      = 1./scl(ind)';
            mom(i).s1 = bsxfun(@times,mom(i).s1,scl1);
            for k=1:K
                mom(i).S2(:,:,k) = (scl1*scl1').*mom(i).S2(:,:,k);
            end
        end
        
        [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);
        
        obj{m}{s}.segment.mom = mom;
        obj{m}{s}.segment.gmm = gmm;
    end
end     
%==========================================================================