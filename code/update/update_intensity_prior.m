function obj = update_intensity_prior(obj,dir_template,iter)
M = numel(obj);
for m=1:M
    pr = do_update(obj{m});    
    
    S = numel(obj{m});
    for s=1:S
        obj{m}{s}.gmm.pr = pr;   
    end
    
    pth1 = fileparts(obj{m}{1}.image(1).fname);
    pth1 = strsplit(pth1,'/');
    pth1 = pth1{end - 1};
    
    fname = fullfile(dir_template,['prior-' pth1 '.mat']);
    save(fname,'pr');
    
    for n=1:size(pr.m,1), fprintf('%2d | pr.m = [%.3f, %s%.3f]\n',iter,pr.m(n,1),sprintf('%.3f, ',pr.m(n,2:end - 1)),pr.m(n,end)); end    
end
%==========================================================================

%==========================================================================
function pr = do_update(obj)
S = numel(obj);

% Compute sample means of all subjects' priors
m0 = 0;
b0 = 0;
n0 = 0;
W0 = 0;
for s=1:S
    m0  = m0 + obj{s}.gmm.pr.m;
    b0  = b0 + obj{s}.gmm.pr.b;
    W0  = W0 + obj{s}.gmm.pr.W;
    n0  = n0 + obj{s}.gmm.pr.n; 
end
m0 = m0/S;
b0 = b0/S;
W0 = W0/S;
n0 = n0/S;

N = size(m0,1);
K = size(m0,2);

for k=1:K
    %______________________________________________________________________________
    %
    % Compute m_0

    g = zeros(N,1);
    H = zeros(N,N);
    for i=1:S    
        [m,b,W,n] = get_po(obj,i);

        g = g + b0(k)*n(k)*W(:,:,k)*(m(:,k)-m0(:,k));
        H = H + b0(k)*n(k)*W(:,:,k);
    end
    m0(:,k) = m0(:,k) + H\g;
    %______________________________________________________________________________


    %______________________________________________________________________________
    %
    % Compute \beta_0

    g_const = 0;
    for i=1:S
        [m,b,W,n] = get_po(obj,i);

        g_const = g_const - 0.5*(N/b(k) + n(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)));
    end
    for subit=1:100,
        % Diff w.r.t. b0
        g  = 0.5*S*N/b0(k) + g_const; % Gradient
        H  = 0.5*S*N/b0(k)^2;         % Hessian
        b0(k) = max(b0(k) + H\g,1e-5);
        if norm(g)==0, break; end
    end
    %______________________________________________________________________________


    %______________________________________________________________________________
    %
    % Set up some constants

    nW_const = zeros(N);
    for i=1:S
        [m,b,W,n] = get_po(obj,i);

        nW_const = nW_const + n(k)*W(:,:,k);
    end

    ElogLam = S*N*log(2);
    for i=1:S
        [m,b,W,n] = get_po(obj,i);

        ElogLam = ElogLam + 2*sum(log(diag(chol(W(:,:,k)))));
        for j=1:N
            ElogLam = ElogLam + psi((n(k)+1-j)/2);
        end
    end

    % convergence = [];
    E = -realmax;

    for it=1:1000
        %______________________________________________________________________________
        %
        % Compute objective function (Equation 10.74 of Bishop)

        oE = E;
        logB = -n0(k)*sum(log(diag(chol(W0(:,:,k))))) - n0(k)*N/2*log(2) - N*(N-1)/4*log(pi);
        for j=1:N, 
            logB = logB - gammaln((n0(k)+1-j)/2); 
        end
        E = (0.5*N*log(b0(k)/(2*pi)) + logB)*S + 0.5*(n0(k)-N-1)*ElogLam;
        for i=1:S
            [m,b,W,n] = get_po(obj,i);

            e = 0.5*(-N*b0(k)/b(k) - b0(k)*n(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)))...
              - 0.5*n(k)*trace(W0(:,:,k)\W(:,:,k));
            E = E + e;
        end
        %if E-oE<abs(E)*eps*D^2, break; end
        if E-oE==0, break; end

    %     convergence = [convergence E];
    %     plot(convergence,'.-'); drawnow;


        %______________________________________________________________________________
        %
        % Compute \nu_0

        % Objective function terms containing n0:
        % NlogB = -n0*N*(sum(log(diag(chol(W0)))) + D/2*log(2));
        % for j=1:D, NlogB = NlogB - N*gammaln((n0+1-j)/2); end
        % E = NlogB + n0*0.5*ElogLam

        g = (sum(log(diag(chol(W0(:,:,k))))) + N/2*log(2))*S - 0.5*ElogLam;
        H = 0;
        for j=1:N
            g = g + S*psi(  (n0(k)+1-j)/2)/2;
            H = H + S*psi(1,(n0(k)+1-j)/2)/4;
        end
        n0(k) = max(n0(k) - H\g,N-0.99999);
        %______________________________________________________________________________



        %______________________________________________________________________________
        %
        % Compute W_0

        % Objective function terms containing W0:
        % E = -n0*N*sum(log(diag(chol(W0))));
        % for i=1:N
        %    E = E - 0.5*n(i)*trace(W0\W(:,:,i));
        % end

        C = inv(chol(W0(:,:,k)));

        % Objective function terms containing W0, after
        % re-expressing using C = inv(chol(W0)):
        % E = n0*N*sum(log(diag(C)));
        % for i=1:N
        %    E = E - 0.5*n(i)*trace(C'*W(:,:,i)*C);
        % end

        G  = -n0(k)*S*diag(1./diag(C)) + nW_const*C;
        for d=1:N,
            c        = C(1:d,d);
            g        = G(1:d,d);
            H        = nW_const(1:d,1:d);
            H(d,d)   = H(d,d) + n0(k)*S/c(d)^2;
            C(1:d,d) = c - H\g;
        end
        C         = inv(C);
        W0(:,:,k) = C'*C;

    end
end

pr.b = b0;
pr.m = m0;
pr.n = n0;
pr.W = W0;
%==========================================================================

%==========================================================================
function [m,b,W,n] = get_po(obj,s)
m = obj{s}.gmm.po.m;
b = obj{s}.gmm.po.b;
W = obj{s}.gmm.po.W;
n = obj{s}.gmm.po.n;
%==========================================================================