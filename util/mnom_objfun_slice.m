function [ll,g,H] = mnom_objfun_slice(r,lnmu,y1,z,lnw)

dm    = size(r);
emu   = cell(dm(4),1);
lnmu1 = cell(dm(4),1);
semu  = 0;
dmu   = cell(dm(4),3);

for k=1:dm(4)
    if nargout>=2
        [lnmu1{k},dmu{k,1},dmu{k,2},dmu{k,3}] = spm_diffeo('bsplins',lnmu(:,:,:,k),y1(:,:,z,:),[1 1 1 0 0 0]);
    else
        lnmu1{k} = spm_diffeo('bsplins',lnmu(:,:,:,k),y1(:,:,z,:),[1 1 1 0 0 0]);
    end
    
    emu{k} = exp(lnmu1{k} + lnw(k));
    semu   = semu + emu{k};
end
for k=1:dm(4)
    emu{k} = emu{k}./semu;
end

if nargout==0, return; end

% Compute gradients and Hessian
msk = isfinite(sum(r(:,:,z,:),4)) & isfinite(semu);
ll  = 0;
if nargout>=2
    g = zeros([dm(1:2),3],'single');
    if nargout>=3
        H = zeros([dm(1:2),6],'single');
    end
end

for k=1:dm(4)
    ftmp  = r(:,:,z,k);
    ll    = ll + sum(sum(ftmp(msk).*log(emu{k}(msk))));

    if nargout>=2
        alpha = emu{k}-r(:,:,z,k);
        for d=1:3
            g(:,:,d) = g(:,:,d) + alpha.*dmu{k,d};
        end

        if nargout>=3
            for k1=1:dm(4)
                if k1~=k,
                    tmp = -emu{k}.*emu{k1};
                else
                    tmp = emu{k}.*(1-emu{k1});
                end

                H(:,:,1)   = H(:,:,1) + tmp.*dmu{k,1}.*dmu{k1,1};
                H(:,:,2)   = H(:,:,2) + tmp.*dmu{k,2}.*dmu{k1,2};
                H(:,:,3)   = H(:,:,3) + tmp.*dmu{k,3}.*dmu{k1,3};
                H(:,:,4)   = H(:,:,4) + tmp.*dmu{k,1}.*dmu{k1,2};
                H(:,:,5)   = H(:,:,5) + tmp.*dmu{k,1}.*dmu{k1,3};
                H(:,:,6)   = H(:,:,6) + tmp.*dmu{k,2}.*dmu{k1,3};
            end
        end
    end
end
if nargout>=2
    g(~isfinite(g)) = 0;
    if nargout>=3
        H(~isfinite(H)) = 0;
    end
end
%==========================================================================