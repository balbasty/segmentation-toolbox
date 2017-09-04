function [g,H,ll] = diff_a(a,r,lnmu,lnw,B,Mmu,Mf,y)
d      = size(r);
Nr     = numel(a);
[E,dE] = spm_dexpm(a,B);
Affine = Mmu\E*Mf;
y1     = affine_transf(Affine,y);
dA     = zeros(4,4,Nr);
for i1=1:Nr
    dA(:,:,i1) = Mmu\dE(:,:,i1)*Mf;
end

lkp = [1 4 5; 4 2 6; 5 6 3];

g = zeros([Nr, 1]);
H = zeros([Nr,Nr]);

ll     = 0;
for z=1:d(3)
    [llz,gz,Hz] = mnom_objfun_slice(r,lnmu,y1,z,lnw);
    ll          = ll + llz;
    
    dAff = cell(Nr,3);
    for i1=1:Nr
        for d1=1:3
            tmp = dA(d1,1,i1)*y(:,:,z,1) + dA(d1,2,i1)*y(:,:,z,2) + dA(d1,3,i1)*y(:,:,z,3) + dA(d1,4,i1);
            dAff{i1,d1} = tmp(:);
        end
    end
    for d1=1:3
        tmp = gz(:,:,d1)';
        tmp = tmp(:)';
        for i1=1:Nr
            g(i1) = g(i1) + double(tmp*dAff{i1,d1});
        end
    end

    % Could probably be speeded up by computing Hessian for changes to elements of Affine,
    % and then transforming this Hessian to obtain changes w.tr.t. rigid parameters.
    for d1=1:3
        for d2=1:3
            tmp1 = Hz(:,:,lkp(d1,d2));
            tmp1 = tmp1(:);
            for i1=1:Nr
                tmp2 = (tmp1.*dAff{i1,d1})';
                for i2=i1:Nr % Only compute upper/lower triangle of symmetric matrix
                    H(i1,i2) = H(i1,i2) + double(tmp2*dAff{i2,d2});
                end
            end
        end
    end
end

% Fill in missing triangle
for i1=1:Nr
    for i2=1:Nr
        H(i2,i1) = H(i1,i2);
    end
end