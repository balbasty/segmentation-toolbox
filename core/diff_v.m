function [g,H,ll] = diff_v(r,lnmu,lnw,Affine,y,int_args)
d   = size(r);
lkp = [1 4 5; 4 2 6; 5 6 3];

g = zeros([d(1:3),3],'single');
H = zeros([d(1:3),6],'single');

y1 = affine_transf(Affine,y);
if int_args > 1
    J  = spm_diffeo('def2jac',y1);
else
    Jz = reshape(Affine(1:3,1:3),[1 1 3 3]);
end

ll  = 0;
for z=1:d(3)
    [llz,gz,Hz] = mnom_objfun_slice(r,lnmu,y1,z,lnw);
    ll          = ll + llz;
    
    if int_args > 1
        Jz = reshape(J(:,:,z,:,:),[d(1:2) 3 3]);
    end

    % Rotate gradients, such that g1 = J'*g0;
    for d1=1:3
        tmp = 0;
        for d2=1:3
            tmp = tmp + Jz(:,:,d2,d1).*gz(:,:,d2);
        end
        g(:,:,z,d1) = tmp;
    end

    % Rotate Hessian, such that H2 = J'*H0*J
    % First do H1 = J'*H0
    RH  = zeros([d(1:2),3,3],'single');
    for d1=1:3
        for d3=1:3
            tmp = 0;
            for d2=1:3
                tmp = tmp + Jz(:,:,d2,d1).*Hz(:,:,lkp(d2,d3));
            end
            RH(:,:,d1,d3) = tmp;
        end
    end

    % Then do H2 = H1*J
    for d1=1:3
        for d3=d1:3 % Only need compute an upper or lower triangle
            tmp = 0;
            for d2=1:3
                tmp = tmp + RH(:,:,d1,d2).*Jz(:,:,d2,d3);
            end
            H(:,:,z,lkp(d1,d3)) = tmp;
        end
    end
end