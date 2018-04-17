function obj = calc_avg_dc(obj)
M = numel(obj);    
for m=1:M
    S      = numel(obj{m});    
    sum_dc = 0;
    for s=1:S    
        % Sum bias field DC components
        sum_dc = sum_dc + obj{m}{s}.segment.bf.dc;
    end
    
    % Average oexp(b1.*b2.*b3.*avg_dc)f bias field DC components
    avg_dc = sum_dc./S;
    
    % Set average bias field DC component
    for s=1:S 
        obj{m}{s}.segment.bf.avg_dc = avg_dc; 
    end
    
    % Convert from log to intensity
    b1 = obj{m}{s}.segment.bf.b1;
    b2 = obj{m}{s}.segment.bf.b2;
    b3 = obj{m}{s}.segment.bf.b3;

    scl = exp(b1.*b2.*b3.*avg_dc);
      
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