function [po,pr,lb] = spm_VBGaussiansFromSuffStats(mom,po,pr)
% Compute estimates of VB-MoG posteriors from sufficient statistics
%
% FORMAT [po,pr,lb] = spm_VBGaussiansFromSuffStats(mom,po,pr)
%
if nargin<2, po = []; end
if nargin<3, pr = []; end

K  = size(mom(1).s0,2);
N  = numel(mom(1).ind);

s0 = 0;
for i=1:numel(mom), 
    s0 = s0 + mom(i).s0; 
end

if isempty(po)
    ib = 0.01*ones(1,K);
    in = N*ones(1,K) - 0.99;
    iW = zeros(N,N,K);
    for k=1:K
        iW(:,:,k) = eye(N);
    end
    im = zeros(N,K);
else
    ib = po.b;
    in = po.n;
    iW = po.W;
    im = po.m;
end

if isempty(pr)
    ib0 = 0.01*ones(1,K);
    in0 = N*ones(1,K) - 0.99;
    iW0 = zeros(N,N,K);
    for k=1:K
        iW0(:,:,k) = eye(N);
    end
    im0 = zeros(N,K);
else
    ib0 = pr.b;
    in0 = pr.n;
    iW0 = pr.W;
    im0 = pr.m;
end

L = -Inf;
for k=1:K,
    b = ib(k);
    n = in(k);
    W = iW(:,:,k);
    m = im(:,k);

    b0 = ib0(k);
    n0 = in0(k);
    W0 = iW0(:,:,k);
    m0 = im0(:,k);

    for iter=1:1024,        
        P  = n*W; 
        C  = inv(P);
        mu = m;        
 
        s1  = zeros(N,1);
        S2  = zeros(N);
        s1i = zeros(N,1);
        S2i = zeros(N,N);

        oL = L;
        ll = 0;
        lq = 0;
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

                % Compute missing data part of objective function
                S  = mom(i).s0(k)*(mux*mux') + mom(i).S2(:,:,k) - mux*mom(i).s1(:,k).' - mom(i).s1(:,k)*mux.';
                ll = ll - 0.5*trace(S/C(ind,ind));
                lq = lq - mom(i).s0(k)*sum(log(diag(chol(C(ind,ind)*2*pi))));  
            end
        end
        

        %----------------------------------------------------------------------
        % Prepare moments, priors and posteriors for further computations
        nmom.s0 = s0(k);
        nmom.s1 = s1;
        nmom.S2 = S2;    

        npr.m = m0;
        npr.b = b0;
        npr.W = W0;
        npr.n = n0;  

        npo.m = m;
        npo.b = b;
        npo.W = W;
        npo.n = n;

        %----------------------------------------------------------------------
        % Compute lower bound
        lb  = lowerbound(nmom,npo,npr);        
        L   = lb + lq + ll;        
%         lb = lb + lq;        
%         fprintf('%d\t%g\n', iter,L);
        if abs((L-oL)/L) < 1e-6, 
            break; 
        end

        if nargout<3
            %----------------------------------------------------------------------
            % Update prior and posterior
            if nargout==2
                [npo,npr] = vmstep(nmom);
                W0        = npr.W;
                m0        = npr.m;
            else
                npo = vmstep(nmom,npr);
            end
           
            m = npo.m;
            n = npo.n;
            W = npo.W;
            b = npo.b;
        end
    end    

    pr.m(:,k)   = m0;
    pr.b(1,k)   = b0;
    pr.W(:,:,k) = W0;
    pr.n(1,k)   = n0;  

    po.m(:,k)   = m;
    po.b(1,k)   = b;
    po.W(:,:,k) = W;
    po.n(1,k)   = n;
end
% fprintf('ML (simple)\t%g\n', ll);
% disp(mu)
% disp(C)