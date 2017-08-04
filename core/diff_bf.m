function [g,H] = diff_bf(chan,bf,r,x,msk,d,K,C,cp,pr,c)

if r==0
    r = ones([d K],'single');
end

r   = reshape(r,[d K]);
msk = reshape(msk,d);
bf  = reshape(bf,d);

d3 = numel(chan(c).T);

H = zeros(d3,d3); 
g = zeros(d3,1); 

for z=1:d(3)

    q = squeeze(r(:,:,z,:));
    q = reshape(q,[prod(d(1:2)) K]);
    
    mskz = msk(:,:,z);
    
    fz = x(:,:,z);    
    fz = fz(mskz);
    
    bfz = bf(:,:,z);
    bfz = bfz(mskz);

    nm = numel(fz);

    cr = cell(C,1);
    for n1=1:C, 
        cr{n1} = double(fz).*double(bfz); 
    end

    w1 = zeros(nm,1);
    w2 = zeros(nm,1);
    for k=1:K
        qk  = q(mskz,k);
        w0  = zeros(nm,1);
        for n1=1:C
            w0 = w0 + pr(n1,c,k)*(cp.po.m(n1,k) - cr{n1});
        end
        w1  = w1 + qk.*w0;
        w2  = w2 + qk*pr(c,c,k);
    end
    wt1       = zeros(d(1:2));
    wt1(mskz) = -(1 + cr{c}.*w1); % US eq. 34 (gradient)
    wt2       = zeros(d(1:2));
    wt2(mskz) = cr{c}.*cr{c}.*w2 + 1; % Simplified Hessian of US eq. 34

    b3 = chan(c).B3(z,:)';
    g  = g + kron(b3,spm_krutil(wt1,chan(c).B1,chan(c).B2,0));
    H  = H + kron(b3*b3',spm_krutil(wt2,chan(c).B1,chan(c).B2,1));
end