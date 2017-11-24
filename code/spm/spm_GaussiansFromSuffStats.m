function varargout = spm_GaussiansFromSuffStats(varargin)
% Compute estimates of Guassian parameters from sufficient statistics
%
% FORMAT varargout = spm_GaussiansFromSuffStats(varargin)
%

% Get/set input
%--------------------------------------------------------------------------
vb     = varargin{1};
mom    = varargin{2};
nomiss = mom(1).nomiss;
if vb
    % Variational Bayes
    %----------------------------------------------------------------------
    if nargin<3
        po = [];
    else
        po = varargin{3};
    end
    
    if nargin<4
        pr = [];
    else
        pr = varargin{4};    
    end        
else
    % Maximum likelihood
    %----------------------------------------------------------------------
    if nargin<3
        vr0 = [];
    else
        vr0 = varargin{3};
    end
    
    if nargin<4
        mn = [];
    else
        mn = varargin{4};    
    end
    
    if nargin<5
        vr = [];
    else
        vr = varargin{5};    
    end
end

K = size(mom(1).s0,2); % Classes
N = numel(mom(1).ind); % Channels

% Sum up 0th moments
s0 = 0;
for i=1:numel(mom), 
    s0 = s0 + mom(i).s0; 
end

% Set initial estimates of parameters
%--------------------------------------------------------------------------
if vb
    % Variational Bayes
    %----------------------------------------------------------------------
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
else
    % Maximum likelihood
    %----------------------------------------------------------------------
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
    
    C1  = 0;
end

% Estimate parameters
%--------------------------------------------------------------------------
ll = zeros(1,K);
lq = zeros(1,K);
lb = zeros(1,K);
for k=1:K,    
    if vb
        % Variational Bayes
        %------------------------------------------------------------------
        b = ib(k);
        n = in(k);
        W = iW(:,:,k);
        m = im(:,k);

        b0 = ib0(k);
        n0 = in0(k);
        W0 = iW0(:,:,k);
        m0 = im0(:,k);
    else
        % Maximum likelihood
        %------------------------------------------------------------------
        C  = iC(:,:,k);
        mu = imu(:,k);
    end
    
    L = -Inf;
    for iter=1:1024,  
        oL = L;
        
        if nomiss            
            i  = numel(mom);
            s1 = mom(i).s1(:,k);
            S2 = mom(i).S2(:,:,k);
        else
            if vb
                % Variational Bayes
                %--------------------------------------------------------------
                mu = m;   
                P  = n*W; 
                C  = inv(P);
            else
                % Maximum likelihood
                %--------------------------------------------------------------
                P  = inv(C);
            end        

            s1    = zeros(N,1);
            S2    = zeros(N);
            s1i   = zeros(N,1);
            S2i   = zeros(N,N);            
            ll(k) = 0;
            lq(k) = 0;
            for i=2:numel(mom),
                if mom(i).s0(k)>eps
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
                    S     = mom(i).s0(k)*(mux*mux') + mom(i).S2(:,:,k) - mux*mom(i).s1(:,k).' - mom(i).s1(:,k)*mux.';
                    ll(k) = ll(k) - 0.5*trace(S/C(ind,ind));
                    lq(k) = lq(k) - mom(i).s0(k)*sum(log(diag(chol(C(ind,ind)*2*pi))));                                                  
    %                 lq(k) = lq(k) - 0.5*mom(i).s0(k)*(LogDet(Pyy) - size(Pyy,1)*(1+log(2*pi)));
                end
            end                    
        end
        
        % Update parameters
        %------------------------------------------------------------------
        if vb
            % Variational Bayes
            %--------------------------------------------------------------
                    
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

            if nargout>=3
                % Update prior and posterior
                if nargout==4
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
            
            lb(k) = lowerbound(nmom,npo,npr);
            L     = lb(k) + lq(k);  
        else
            % Maximum likelihood
            %--------------------------------------------------------------                                    
            if nargout>=3
                mu = s1/s0(k);        
                C  = (S2 - s1*s1'/s0(k) + N*vr0)/(s0(k) + N);        
            end
            
%             ll(k) = -0.5*s0(k)*(LogDet(C) + N*log(2*pi)) - 0.5*s0(k)*trace((mu*mu' + S2/s0(k) - 2*mu*s1'/s0(k))/C);
            L = ll(k) + lq(k);
        end
        
        % Check convergence
        %--------------------------------------------------------------------------           
%         fprintf('%d\t%g\n', iter,L);
        if abs((L - oL)/L) < 1e-6 || nomiss, 
            % Finished
            break; 
        end
    end    

    if vb
        % Variational Bayes
        %------------------------------------------------------------------
        pr.m(:,k)   = m0;
        pr.b(1,k)   = b0;
        pr.W(:,:,k) = W0;
        pr.n(1,k)   = n0;  

        po.m(:,k)   = m;
        po.b(1,k)   = b;
        po.W(:,:,k) = W;
        po.n(1,k)   = n;
    else
        % Maximum likelihood
        %------------------------------------------------------------------
        mn(:,k)     = mu;
        vr(:,:,k)   = C;
        C1          = C1 + (S2 - s1*s1'/s0(k));
    end
end

% Prepare output
%--------------------------------------------------------------------------
varargout{1} = sum(lb) + sum(lq);
if vb     
    % Variational Bayes
    %----------------------------------------------------------------------        
    varargout{2} = s0;
    varargout{3} = po;
    varargout{4} = pr;
else
    % Maximum likelihood
    %----------------------------------------------------------------------
    vr1 = zeros(N,N,K);
    for k=1:K,
        vr1(:,:,k) = (C1 + N*vr0)/(sum(s0) + N);
    end 
    
    varargout{2} = s0;
    varargout{3} = mn;
    varargout{4} = vr;
    varargout{5} = vr1;
end