function varargout = update_bf(varargin)
% Update bias field parameters (in chan{:}.T)
% FORMAT varargout = update_bf(varargin)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging


% Read function input
%--------------------------------------------------------------------------
ll       = varargin{1};
llrb     = varargin{2};
llr      = varargin{3};
buf      = varargin{4};
mg       = varargin{5};
gmm      = varargin{6};
wp       = varargin{7};
lkp      = varargin{8};
chan     = varargin{9};
fig      = varargin{10};
L        = varargin{11};
print_ll = varargin{12};
armijo   = varargin{13};
wp_l     = varargin{14};

% Some parameters
nz = numel(buf);
N  = numel(buf(1).f);
K  = numel(mg);
d  = size(buf(1).msk{1});

% Get means and precisions
pr = zeros(size(gmm.po.W));
for k=1:K 
    pr(:,:,k) = gmm.po.n(k)*gmm.po.W(:,:,k); 
end
mn = gmm.po.m; 

for n=1:N
    d3  = numel(chan(n).T);
    if d3>0
        % Compute objective function and its first and second derivatives
        Alpha = zeros(d3,d3); % Second derivatives
        Beta  = zeros(d3,1);  % First derivatives

        for z=1:nz
            if ~buf(z).nm(n), continue; end

            cr                 = cell(N,1);
            for n1=1:N, cr{n1} = double(buf(z).f{n1}).*double(buf(z).bf{n1}); end

            q = latent(buf(z).f,buf(z).bf,mg,gmm,buf(z).dat,lkp,wp,buf(z).msk,buf(z).code,buf(z).labels,wp_l,cr);

            w1 = zeros(buf(z).nm(n),1);
            w2 = zeros(buf(z).nm(n),1);
            for k=1:K
                qk  = q(buf(z).msk{n},k);
                w0  = zeros(prod(d(1:2)),1);
                for n1=1:N
                    w0(buf(z).msk{n1}) = w0(buf(z).msk{n1}) + pr(n1,n,k)*(mn(n1,k) - cr{n1});
                end
                w1  = w1 + qk.*w0(buf(z).msk{n});
                w2  = w2 + qk*pr(n,n,k);
            end
            wt1                = zeros(d(1:2));
            wt1(buf(z).msk{n}) = -(1 + cr{n}.*w1); % US eq. 34 (gradient)
            wt2                = zeros(d(1:2));
            wt2(buf(z).msk{n}) = cr{n}.*cr{n}.*w2 + 1; % Simplified Hessian of US eq. 34
            clear cr

            b3    = chan(n).B3(z,:)';
            Beta  = Beta  + kron(b3,spm_krutil(wt1,chan(n).B1,chan(n).B2,0));
            Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,chan(n).B1,chan(n).B2,1));
            clear wt1 wWritet2 b3
        end

        oll     = ll;
        C       = chan(n).C; % Inverse covariance of priors
        oldT    = chan(n).T;

        % Gauss-Newton update of bias field parameters
        Update  = reshape((Alpha + C)\(Beta + C*chan(n).T(:)),size(chan(n).T));
        clear Alpha Beta
        
        for line_search=1:12
            chan(n).T = chan(n).T - armijo*Update; % Backtrack if necessary

            % Re-generate bias field, and compute terms of the objective function
            chan(n).ll = double(-0.5*chan(n).T(:)'*C*chan(n).T(:));
            for z=1:nz
                if ~buf(z).nm(n), continue; end
                
                bf           = transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T);
                tmp          = bf(buf(z).msk{n});
                buf(z).bf{n} = single(exp(tmp));
            end

            ollrb            = llrb;
            llrb             = 0;
            for n1=1:N, llrb = llrb + chan(n1).ll; end
            ll               = llr + llrb;

            % Compute responsibilities and moments
            [mom,dll] = compute_moments(buf,lkp,mg,gmm,wp,wp_l);        
            ll        = ll + dll; 

            % Compute missing data and VB components of ll
            dll  = spm_VBGaussiansFromSuffStats(mom,gmm);
            ll   = ll + sum(sum(dll));

            if ll>=oll
                L{1}(end + 1) = ll;
                L{2}(end + 1) = llrb;
                armijo        = min(armijo*1.25,1);
                debug_view('convergence',fig{4},lkp,buf,L);
                my_fprintf('Bias-%d:\t%g\t%g\t%g :o)\n', n, ll, llr,llrb,print_ll);
                break;
            else
                ll        = oll;
                llrb      = ollrb;
                chan(n).T = oldT;
                armijo    = armijo*0.5;
                my_fprintf('Bias-%d:\t%g\t%g\t%g :o(\n', n, ll, llr,llrb,print_ll);
                if line_search==12
                    L{1}(end + 1) = ll;
                    L{2}(end + 1) = llrb;
                end                    
            end 
        end
        clear oldT
    end
end

% Write function output
%--------------------------------------------------------------------------
varargout{1} = ll;
varargout{2} = llrb;
varargout{3} = buf;
varargout{4} = chan;
varargout{5} = L;
varargout{6} = armijo;
%========================================================================== 

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%==========================================================================