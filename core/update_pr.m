function pr = update_pr(vbmog)

itout = 50;
itin  = 20;

pr = vbmog(1).pr;
po = {};
for n=1:numel(vbmog)
    po{n} = vbmog(n).po;
end

In.C = size(pr.m,1);
In.K = size(pr.m,2);

s = 1:numel(po);

M = numel(s);

for t=1:In.K
    
    lLambda = 0;
    
    for subj=s(1):s(end)
        lLambda = lLambda + Elogdet(po{subj}.W(:,:,t),po{subj}.nu(t));
    end
    
    for it=1:itout
        
        %terms of objective function that don't depend on posteriors
        E   = 0.5*M*(In.C*log(pr.beta(t)/(2*pi)) + 2*BWishartLog(pr.W(:,:,t),pr.nu(t)))  + 0.5*(pr.nu(t) - In.C)*lLambda;
        
        %gradient with respect to pr.beta
        g_b = 0;
        
        %Hessian with respect to pr.m
        H_m = zeros(In.C,In.C);
        
        %part of derivatives with respect to pr.W
        nW = zeros(In.C);
        for subj=s(1):s(end)
            nW = nW + po{subj}.nu(t)*po{subj}.W(:,:,t);
        end
        
        %part of objective function (depending on subject posteriors)
        for subj=s(1):s(end)
            d       = po{subj}.m(:,t) - pr.m(:,t);
            E       = E + 0.5*(- In.C*pr.beta(t)/po{subj}.beta(t) - po{subj}.nu(t)*...
                trace(pr.W(:,:,t)\po{subj}.W(:,:,t) + pr.beta(t)*(d*d')*po{subj}.W(:,:,t)));
            g_b     = g_b - 0.5*(In.C/po{subj}.beta(t) + po{subj}.nu(t)*d'*po{subj}.W(:,:,t)*d);
        end
        
        %update pr.beta
        for subit=1:itin
            g           = 0.5*M*In.C/pr.beta(t) + g_b;
            H           = -0.5*M*In.C/pr.beta(t)^2;
            pr.beta(t) = max(pr.beta(t) - H\g,1e-3);
            if norm(g)==0, break; end
        end
        
        %update pr.m
        for subj=s(1):s(end)
            H_m = H_m + pr.beta(t)*po{subj}.nu(t)*po{subj}.W(:,:,t);
        end
        for subit=1:itin
            g = zeros(In.C,1);
            for subj=s(1):s(end)
                g   = g + pr.beta(t)*po{subj}.nu(t)*po{subj}.W(:,:,t)*(po{subj}.m(:,t)-pr.m(:,t));
            end
            pr.m(:,t) = max(pr.m(:,t) + (H_m\g),1e-3);
        end
        
        %update pr.nu
        for subit=1:itin
            g = (sum(log(diag(chol(pr.W(:,:,t))))) + In.C/2*log(2))*M - 0.5*lLambda;
            H = 0;
            for j=1:In.C
                g = g + M*psi(  (pr.nu(t)+1-j)/2)/2;
                H = H + M*psi(1,(pr.nu(t)+1-j)/2)/4;
            end
            pr.nu(t) = max(pr.nu(t) - H\g,In.C-0.99999);
        end
        
        %update pr.W
        for subit=1:itin
            Ck = inv(chol(pr.W(:,:,t)));
            G  = -pr.nu(t)*M*diag(1./diag(Ck)) + nW*Ck;
            for d=1:In.C,
                c         = Ck(1:d,d);
                g         = G(1:d,d);
                H         = nW(1:d,1:d);
                H(d,d)    = H(d,d) + pr.nu(t)*M/c(d)^2;
                Ck(1:d,d) = c - H\g;
            end
            Ck  = inv(Ck);
            pr.W(:,:,t) = Ck'*Ck;
        end
    end
end

%==========================================================================
function logB=BWishartLog(W,nu)
%Logarithm of normalization constant for Wishart Distribution

D         = size(W,1);
prodgamln = sum(gammaln(0.5*(nu + 1 - [1:D])));
logB      = - 0.5*nu*logdet(W) - 0.5*nu*D*log(2) - 0.25*D*(D-1)*log(pi) - prodgamln;
%==========================================================================

%==========================================================================
function lLambda = Elogdet(W,nu)
%Expectation of logarithm of determinant of Wishart distributed covariance matrix

C       = size(W,1);
ldW     = 2*sum(log(diag(chol(W))));
Sdig    = sum(psi((1 + nu-(1:C))/2));
lLambda = ldW + C*log(2) + Sdig;
%==========================================================================

%==========================================================================
function l = logdet(W)

l = 2*sum(log(diag(chol(W))));
%==========================================================================
