function [mg,mn,vr] = spm_GaussiansFromSuffStats(mom)
% Compute ML estimates of Gaussians from sufficient statistics
%
% FORMAT [mg,mn,vr] = spm_GaussiansFromSuffStats(mom)
%

K  = size(mom(1).s0,2);
N  = numel(mom(1).ind);
mg = zeros(1,K);
mn = zeros(N,K);
vr = zeros(N,N,K);
ll = -Inf;
for k=1:K,
    C  = eye(N);
    mu = zeros(N,1);
    s0 = 0;
    for i=1:numel(mom), s0 = s0 + mom(i).s0; end
    for iter=1:1024,
        s1     = zeros(N,1);
        S2     = zeros(N);
        P      = inv(C);
        s1i    = zeros(N,1);
        S2i    = zeros(N,N);
        old_ll = ll;
        ll     = 0;
        for i=1:numel(mom),
            if mom(i).s0,
                ind            = mom(i).ind;
                mux            = mu( ind);
                muy            = mu(~ind);
                Pyy            = P(~ind,~ind);
                Pyx            = P(~ind, ind);
                R              = Pyy\Pyx;
                Ex             = mom(i).s1/mom(i).s0;
                Exx            = mom(i).S2/mom(i).s0;
                tmp            = R*(mux-Ex)*muy';
                s1i( ind)      = mom(i).s1;
                S2i( ind, ind) = mom(i).S2;
                s1i(~ind)      = mom(i).s0*(R*(mux-Ex) + muy);
                S2i(~ind,~ind) = mom(i).s0*(R*(mux*mux'-mux*Ex'-Ex*mux'+Exx)*R' + tmp + tmp' + muy*muy' + inv(Pyy));
                S2i( ind,~ind) = mom(i).s0*(R*(mux*Ex'-Exx) + muy*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1             = s1 + s1i;
                S2             = S2 + S2i;

                % Compute objective function
                S              = mom(i).s0*(mux*mux') + mom(i).S2 - mux*mom(i).s1.' - mom(i).s1*mux.';
                ll             = ll - mom(i).s0*sum(log(diag(chol(C(ind,ind)*2*pi)))) - 0.5*trace(S/C(ind,ind));
            end
        end
        mu = s1/s0;
        C  = (S2 - s1*s1'/s0)/(s0+eps);
        %fprintf('%d\t%g\n', iter,ll);
        if ll-old_ll < 1e-12, break; end
    end
    mg(1,k)   = s0;
    mn(:,k)   = mu;
    vr(:,:,k) = C;
end
fprintf('ML (simple)\t%g\n', ll);
disp(mu)
disp(C)

