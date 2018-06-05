function [mom,ll,mrf,mgm] = compute_moments(buf,lkp,mg,gmm,wp,wp_l,resp_read,resp_save,mrf,update_ElnPzN)
if nargin<10, update_ElnPzN = true; end
    
K   = numel(lkp.part);
Kb  = max(lkp.part);
nz  = numel(buf);
N   = numel(buf(1).f);
dm  = size(buf(1).msk{1});
if nz==1, dm = [dm 1]; end
I   = prod(dm(1:2));
mgm = zeros(1,Kb);

if mrf.do_mrf        
    % Compute neighborhood part of responsibilities
    % (from previous responsibilities)
    lnPzN = compute_lnPzN(mrf,lkp,resp_read);
end

% Compute responsibilities
ll  = 0;
mom = moments_struct(K,N);
for z=1:nz
    if ~buf(z).Nm, continue; end
        
    if nargout==4
        % For updating tissue weights
        s   = 1./(double(buf(z).dat)*wp');
        mgm = mgm + s'*double(buf(z).dat);
        clear s
    end
    
    % Get BX (bias-field x image)
    cr                             = NaN(I,N);
    for n=1:N, cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end    
       
    % Get neighborhood term
    if mrf.do_mrf   
        lnPzN1 = double(reshape(lnPzN(:,:,z,:),[I Kb]));
    else
        lnPzN1 = zeros([1 Kb]);
    end
    
    % Compute VB-responsibilities
    [q,dll] = latent(buf(z).f,buf(z).bf,mg,gmm,double(buf(z).dat),lnPzN1,lkp,wp,buf(z).msk,buf(z).code,buf(z).labels,wp_l,cr);
    ll      = ll + dll;
    clear lnPzN1
    
    % Store responsibilities in a NIfTI file
    for k=1:K
        tmp                     = single(reshape(sum(q(:,k),2),dm));
        tmp(~isfinite(tmp))     = 0;
        resp_save(k).dat(:,:,z) = tmp;
    end  
    clear tmp   
    
    % Update sufficient statistics
    mom = spm_SuffStats(cr,q,mom,buf(z).code);
end

if mrf.do_mrf && update_ElnPzN       
    % Compute neighborhood part of lower bound 
    % (always requires most updated responsibilities)
    mrf.ElnPzN  = compute_ElnPzN(mrf,lkp,resp_save);    
    ll          = ll + mrf.ElnPzN;
else
    ll          = ll + mrf.ElnPzN;
end
%==========================================================================

%==========================================================================
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

%==========================================================================
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

lnPzN  = spm_vbmrf_lowerbound(R,single(mrf.ElnUpsilon),mrf.w);
ElnPzN = sum(sum(sum(sum(lnPzN))));
%==========================================================================