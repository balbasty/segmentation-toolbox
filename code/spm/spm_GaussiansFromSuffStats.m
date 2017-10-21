function [mg,mn,vr] = spm_GaussiansFromSuffStats(mom,vr0)
% Compute ML estimates of Gaussians from sufficient statistics
%
% FORMAT [mg,mn,vr] = spm_GaussiansFromSuffStats(mom)
%
if nargin<2, vr0 = 0; end

tiny = eps*eps;

K  = size(mom(1).s0,2);
N  = numel(mom(1).ind);
mg = zeros(1,K);
mn = zeros(N,K);
vr = zeros(N,N,K);
ll = -Inf;

s0 = zeros(K,1);
for k=1:K
    for i=1:numel(mom), 
        s0(k) = s0(k) + mom(i).s0(k); 
    end
end

for k=1:K,
    C  = eye(N);
    mu = zeros(N,1);    
    for iter=1:1024,
        s1     = zeros(N,1);
        S2     = zeros(N);
        P      = inv(C);
        s1i    = zeros(N,1);
        S2i    = zeros(N,N);
        old_ll = ll;
        ll     = 0;
%         llz    = 0;
        for i=1:numel(mom),
            if mom(i).s0(k),
                ind            = mom(i).ind;
                mux            = mu( ind,:);
                muy            = mu(~ind,:);
                Pyy            = P(~ind,~ind);
                Pyx            = P(~ind, ind);
                R              = Pyy\Pyx;
                Ex             = mom(i).s1(:,k)/mom(i).s0(k);
                Exx            = mom(i).S2(:,:,k)/mom(i).s0(k);
                tmp            = R*(mux-Ex)*muy';
                s1i( ind)      = mom(i).s1(:,k);
                S2i( ind, ind) = mom(i).S2(:,:,k);
                s1i(~ind)      = mom(i).s0(k)*(R*(mux-Ex) + muy);
                S2i(~ind,~ind) = mom(i).s0(k)*(R*(mux*mux'-mux*Ex'-Ex*mux'+Exx)*R' + tmp + tmp' + muy*muy' + inv(Pyy));
                S2i( ind,~ind) = mom(i).s0(k)*(R*(mux*Ex'-Exx) + muy*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1             = s1 + s1i;
                S2             = S2 + S2i;

                % Compute objective function
                S   = mom(i).s0(k)*(mux*mux') + mom(i).S2(:,:,k) - mux*mom(i).s1(:,k).' - mom(i).s1(:,k)*mux.';
                ll  = ll - mom(i).s0(k)*sum(log(diag(chol(C(ind,ind)*2*pi)))) - 0.5*trace(S/C(ind,ind));                
%                 llz = llz + 0.5*mom(i).s0(k)*(LogDet(Pyy) - size(Pyy,1)*(1+log(2*pi)));
            end
        end
        mu = s1/(s0(k) + tiny);
        C  = (S2 - s1*s1'/s0(k) + N*vr0)/(s0(k) + N);
        %fprintf('%d\t%g\n', iter,ll);
        if ll-old_ll < 1e-12, break; end
    end
    mg(k)     = s0(k);
    mn(:,k)   = mu;
    vr(:,:,k) = C;
end
% fprintf('ML (simple)\t%g\n', ll);
% disp(mu)
% disp(C)
%==========================================================================