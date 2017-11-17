function [mn,vr,vr1] = spm_GaussiansFromSuffStats(mom,vr0,mn,vr)
% Compute ML estimates of Gaussians from sufficient statistics
%
% FORMAT [mn,vr,vr1] = spm_GaussiansFromSuffStats(mom)
%
if nargin<3, mn = []; end
if nargin<4, vr = []; end

K = size(mom(1).s0,2);
N = numel(mom(1).ind);

s0 = 0;
for i=1:numel(mom), 
    s0 = s0 + mom(i).s0; 
end

if isempty(mn)
    imu = zeros(N,K);
    mn  = zeros(N,K);
else
    imu = mn;
end

if isempty(vr)
    iC = zeros(N,N,K);
    for k=1:K
        iC(:,:,k) = eye(N);
    end
    vr = zeros(N,N,K);
else
    iC = vr;
end

vr1 = zeros(N,N,K);
C1  = 0;
L   = -Inf;
for k=1:K
    C  = iC(:,:,k);
    mu = imu(:,k);
    for iter=1:1024
        s1  = zeros(N,1);
        S2  = zeros(N);
        P   = inv(C);
        s1i = zeros(N,1);
        S2i = zeros(N,N);
        oL = L;
        ll = 0;
        lq = 0;
        for i=1:numel(mom)
            if mom(i).s0(k)
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
                S  = mom(i).s0(k)*(mux*mux') + mom(i).S2(:,:,k) - mux*mom(i).s1(:,k).' - mom(i).s1(:,k)*mux.';
                ll = ll - 0.5*trace(S/C(ind,ind));
                lq = lq - mom(i).s0(k)*sum(log(diag(chol(C(ind,ind)*2*pi))));                
            end
        end
        mu = s1/s0(k);        
        C  = (S2 - s1*s1'/s0(k) + N*vr0)/(s0(k) + N);        
        
        L  = ll + lq;        
%         fprintf('%d\t%g\n', iter,L);
        if abs((L-oL)/L) < 1e-6, 
            break; 
        end
    end
    mn(:,k)   = mu;
    vr(:,:,k) = C;
    C1        = C1 + (S2 - s1*s1'/s0(k));
end
for k=1:K,
    vr1(:,:,k) = (C1 + N*vr0)/(sum(s0) + N);
end
% fprintf('ML (simple)\t%g\n', ll);
% disp(mu)
% disp(C)