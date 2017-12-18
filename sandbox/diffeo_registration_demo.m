function diffeo_registration_demo

addpath(genpath('../code'));

Nmu = nifti('/home/mbrud/Dropbox/PhD/Data/JA/diffeo_registration_demo/Template_4.nii');
Nf  = nifti(char('/home/mbrud/Dropbox/PhD/Data/JA/diffeo_registration_demo/rc1mprage0040063.nii','/home/mbrud/Dropbox/PhD/Data/JA/diffeo_registration_demo/rc2mprage0040063.nii'));

mu  = single(log(Nmu.dat(:,:,:,:)));
f   = single(cat(4,Nf(1).dat(:,:,:),Nf(2).dat(:,:,:)));
f   = cat(4,f,max(1-sum(f,4),0));

dm       = size(f);
w        = zeros(dm(4),1); % This should contain logs of mixing proportions
rparam   = [0 0.005 1 0.25 1]; % Decrease regularisation over iterations 
int_args = 8;                  % Set to 1 for small deformation approximation

% Basis functions for the lie algebra of the special Eucliden group
% (SE(3)).

B = affine_basis;

Nr       = size(B,3);
y        = cell(3,1);
[y{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
y        = cat(4,y{:});
vx       = sqrt(sum(Nf(1).mat(1:3,1:3).^2));
prm      = [vx rparam];
if int_args > 1
    Greens = spm_shoot_greens('kernel',dm(1:3),prm);
else
    Greens = [];
end

v        = zeros([dm(1:3),3],'single');
r        = zeros([Nr,1],'single');

E      = spm_dexpm(r,B);
Affine = Nmu.mat\E*Nf(1).mat;
y1     = affine_transf(Affine,y);
ll     = [mnom_objfun(f,mu,y1,w) 0];
fprintf('%d %g %g %g\n', 0, -ll, -sum(ll));

for it=1:30
    oll            = ll;
    [r,ll1,Affine] = update_affine(r,f,mu,w,B,Nmu.mat,Nf(1).mat,y);
    ll(1)          = ll1;
    fprintf('%d %g %g %g\n', it, -ll, -sum(ll));
    [v,ll,y]       = update_nonlin(v,f,mu,w,Affine,y,prm,int_args,Greens);
    fprintf('%d %g %g %g\n', it, -ll, -sum(ll));
    if abs(sum(ll-oll)/sum(ll))<1e-3
        break;
    end
end



function [v,ll,y] = update_nonlin(v,f,mu,w,Affine,y,prm,int_args,Greens)
dm     = size(f);
lkp    = [1 4 5; 4 2 6; 5 6 3];

g   = zeros([dm(1:3),3],'single');
H   = zeros([dm(1:3),6],'single');
y1  = affine_transf(Affine,y);
if int_args > 1
    J  = spm_diffeo('def2jac',y1);
else
    Jz = reshape(Affine(1:3,1:3),[1 1 3 3]);
end
u   = spm_diffeo('vel2mom',v,prm);
ll  = [0 -0.5*sum(sum(sum(sum(u.*v))))];
for z=1:dm(3)
    [llz,gz,Hz] = mnom_objfun_slice(f,mu,y1,z,w);
    ll(1)       = ll(1) + llz;

    if int_args > 1
        Jz = reshape(J(:,:,z,:,:),[dm(1:2) 3 3]);
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
    RH  = zeros([dm(1:2),3,3],'single');
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

g  = g + u;
dv = spm_diffeo('fmg',H, g, [prm 3 2]);

if false
    v   = v - dv;
    if nargout>=3
        y     = make_deformation(v,prm,int_args, Greens);
    end
else
    scale   = 1;
    oll     = ll;
    ov      = v;
    for subit=1:4
        v     = v - scale*dv;
        y     = make_deformation(v,prm,int_args, Greens);
        y1    = affine_transf(Affine,y);
        ll(1) = mnom_objfun(f,mu,y1,w);
        u     = spm_diffeo('vel2mom',v,prm);
        ll(2) = -0.5*sum(sum(sum(sum(u.*v))));
        if sum(ll)>sum(oll)
           %scale = min(scale*1.25,1);
            break;
        else
            scale = scale*0.5;
            v     = ov;
            y     = make_deformation(v,prm,int_args, Greens);
        end
    end

    if sum(ll)<sum(oll)
        ll = oll;
        v  = ov;
    end
end


function y = make_deformation(v,prm,int_args, K)
dm        = size(v);
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
id        = cat(4,id{:});
if int_args<=1
    y = id + v;
else
    y   = spm_shoot3d(v,prm,int_args, K);
end



function [r,ll,Affine] = update_affine(r,f,mu,w,B,Mmu,Mf,y)
dm     = size(f);
Nr     = numel(r);
[E,dE] = spm_dexpm(r,B);
Affine = Mmu\E*Mf;
y1     = affine_transf(Affine,y);
dA     = zeros(4,4,Nr);
for i1=1:Nr
    dA(:,:,i1) = Mmu\dE(:,:,i1)*Mf;
end
lkp    = [1 4 5; 4 2 6; 5 6 3];
scale  = 1;
ga     = zeros([Nr, 1],'single');
Ha     = zeros([Nr,Nr],'single');
ll     = 0;
for z=1:dm(3)
    [llz,gz,Hz] = mnom_objfun_slice(f,mu,y1,z,w);
    ll          = ll + llz;

    dAff        = cell(Nr,3);
    for i1=1:Nr
        for d1=1:3
            tmp = dA(d1,1,i1)*y(:,:,z,1) + dA(d1,2,i1)*y(:,:,z,2) + dA(d1,3,i1)*y(:,:,z,3) + dA(d1,4,i1);
            dAff{i1,d1} = tmp(:);
        end
    end
    for d1=1:3
        tmp = gz(:,:,d1);
        tmp = tmp(:)';
        for i1=1:Nr
            ga(i1) = ga(i1) + tmp*dAff{i1,d1};
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
                    Ha(i1,i2) = Ha(i1,i2) + tmp2*dAff{i2,d2};
                end
            end
        end
    end
end

% Fill in missing triangle
for i1=1:Nr
    for i2=1:Nr
        Ha(i2,i1) = Ha(i1,i2);
    end
end
dr      = Ha\ga;

if false
    r      = r - dr;
    E      = spm_dexpm(r,B);
    Affine = Mmu\E*Mf;
else
    % Linesearch
    oll     = ll;
    or      = r;
    for subit=1:4
        r      = r - scale*dr;
        E      = spm_dexpm(r,B);
        Affine = Mmu\E*Mf;
        y1     = affine_transf(Affine,y);
        ll1    = mnom_objfun(f,mu,y1,w);
        ll     = ll1;
        if ll>oll
           %scale = min(scale*1.2,1);
            break;
        else
            scale  = scale*0.5;
            r      = or;
        end
    end

    if ll<oll
        ll = oll;
        r  = or;
    end
end



function ll = mnom_objfun(f,mu,y1,w)
dm   = size(f);
ll   = 0;
for z=1:dm(3)
    ll = ll + mnom_objfun_slice(f,mu,y1,z,w);
end



function [ll,g,H] = mnom_objfun_slice(f,mu,y1,z,w)
dm   = size(f);
emu  = cell(dm(4),1);
mu1  = cell(dm(4),1);
semu = 0;
dmu  = cell(dm(4),3);
for k=1:dm(4)
    if nargout>=2
        [mu1{k},dmu{k,1},dmu{k,2},dmu{k,3}] = spm_diffeo('bsplins',mu(:,:,:,k),y1(:,:,z,:),[1 1 1  0 0 0]);
    else
        mu1{k} = spm_diffeo('bsplins',mu(:,:,:,k),y1(:,:,z,:),[1 1 1  0 0 0]);
    end
    emu{k} = exp(mu1{k} + w(k));
    semu   = semu + emu{k};
end
for k=1:dm(4)
    emu{k} = emu{k}./semu;
    if z ==round((dm(3)+1)/2);
        subplot(2,dm(4),k      ); imagesc(emu{k}');     axis image xy off
        subplot(2,dm(4),k+dm(4)); imagesc(f(:,:,z,k)'); axis image xy off
    end
end
drawnow

if nargout==0, return; end

msk = isfinite(f(:,:,z,1)) & isfinite(semu);
ll  = 0;
if nargout>=2
    g   = zeros([dm(1:2),3],'single');
    if nargout>=3
        H   = zeros([dm(1:2),6],'single');
    end
end

for k=1:dm(4)
    ftmp  = f(:,:,z,k);
    ll    = ll + sum(sum(ftmp(msk).*log(emu{k}(msk))));

    if nargout>=2
        alpha = emu{k}-f(:,:,z,k);
        for d=1:dm(4)
            g(:,:,d) = g(:,:,d) + alpha.*dmu{k,d};
        end

        if nargout>=3
            for k1=1:dm(4)
                if k1~=k,
                    tmp = -emu{k}.*emu{k1};
                else
                    tmp = max(emu{k}.*(1-emu{k1}),0);
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


function y1 = affine_transf(Affine,y)
y1     = cat(4, Affine(1,1)*y(:,:,:,1) + Affine(1,2)*y(:,:,:,2) + Affine(1,3)*y(:,:,:,3) + Affine(1,4),...
                Affine(2,1)*y(:,:,:,1) + Affine(2,2)*y(:,:,:,2) + Affine(2,3)*y(:,:,:,3) + Affine(2,4),...
                Affine(3,1)*y(:,:,:,1) + Affine(3,2)*y(:,:,:,2) + Affine(3,3)*y(:,:,:,3) + Affine(3,4));

