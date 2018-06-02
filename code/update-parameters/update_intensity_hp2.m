function obj = update_intensity_hp2(obj)
% FORMAT obj = update_intensity_hp2(obj)
% 
% Hierarchical Gauss-Wishart intensity prior, with a Wishart hyper-prior on
% prior Precision matrices to make them similar.
% It should bias the GMM segmentation towards a "K-means" behaviour (were
% all variances are the same).

% M = number of modality families
% N = number of channels / family
% S = number of subjects
% K = number of Gaussians in the mixture

dir_template = obj{1}{1}.dir_template;
M            = numel(obj);

obj_ct = {};
for m=1:M
    if strcmp(obj{m}{1}.modality,'CT')
        obj_ct{end + 1} = obj{m};        
    end
end

if ~isempty(obj_ct)
    % CT update
    %----------------------------------------------------------------------
    
    obj_tmp = {};
    for m=1:numel(obj_ct)
        S = numel(obj_ct{m});                
        for s=1:S
            if obj_ct{m}{s}.status==0
                obj_tmp{end + 1}.gmm = obj_ct{m}{s}.segment.gmm;
            end
        end
    end
    
    pr = do_update(obj_tmp,obj_ct{1}{1}.segment.constr_inthp);

    for m=1:M
        if strcmp(obj{m}{1}.modality,'CT')
            S = numel(obj{m}); 
            for s=1:S            
                obj{m}{s}.segment.gmm.pr = pr;               
            end
        end
    end    

    fname = fullfile(dir_template,['prior-ct.mat']);
    save(fname,'pr');
    
    clear obj_tmp
end

% MRI update
%--------------------------------------------------------------------------
for m=1:M
    if strcmp(obj{m}{1}.modality,'MRI')
        S       = numel(obj{m});
        obj_tmp = {};
        for s=1:S
            if obj{m}{s}.status==0
                obj_tmp{end + 1}.gmm = obj{m}{s}.segment.gmm;
            end
        end

        pr = do_update(obj_tmp,obj{m}{1}.segment.constr_inthp);    

        for s=1:S
            obj{m}{s}.segment.gmm.pr = pr;   
        end

        fname = fullfile(dir_template,['prior-' num2str(m) '.mat']);
        save(fname,'pr');
    end
end
%==========================================================================


%==========================================================================
function pr = do_update(obj,constrained,tol)
% FORMAT pr = do_update(obj,consrained,tol)
%
% obj - cell of size S (the number of subjects) containing structures with
%       the field gmm.pr
%       > pr is a structure containing posterior parameters of the
%         Gauss-Wishart distribution
%
% constrained - If true:  constrain covariances to be alike
%               If false: mode estimate for all parameters
%
% tol - stopping tolerance [1e-4]
%
% pr  - structure with fields
%       * m0, b0, W0, n0 (Gauss-Wishart prior)     -> length K
%       if constrained:
%       * p, V           (hyper-Wishart posterior) -> length K
%       * p0, V0         (hyper-Wishart prior)     -> length 1
%
% Hierarchical Gauss-Wishart intensity prior, with a Wishart hyper-prior on
% all prior scale matrices to make them similar.
% It should bias the GMM segmentation towards a "K-means" behaviour (were
% all variances are the same).

% M = number of modality families
% N = number of channels / family
% S = number of subjects
% K = number of Gaussians in the mixture
if nargin<3, tol = 1e-4; end

S = numel(obj);

% -------------------------------------------------------------------------
% Starting estimates
% > Compute sample means of all posterior Gauss-Wishart parameters
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

% preallocate
LogDetW0  = zeros(size(n0));
V         = zeros(size(W0));
p         = zeros(size(n0));
p0        = 0;

part = obj{s}.gmm.pr.part;

% -------------------------------------------------------------------------
%   Gauss-Wishart "mean" parameters
% -------------------------------------------------------------------------

for k=1:K
    
    % ---------------------------------------------------------------------
    % Update m0 (mode, closed-form)
    
    Lambda   = 0;
    LambdaMu = 0;
    for s=1:S
        [m,~,W,n] = get_po(obj,s);
        Lambda    = Lambda   + n(k)*W(:,:,k);
        LambdaMu  = LambdaMu + n(k)*W(:,:,k)*m(:,k);
    end
    m0(:,k) = Lambda \ LambdaMu;
    
    % ---------------------------------------------------------------------

    
    % ---------------------------------------------------------------------
    % Update b0 (mode, closed-form)

    b0(k)= 0;
    for s=1:S
        [m,b,W,n] = get_po(obj,s);
        m1 = m(:,k) - m0(:,k);
        b0(k) = b0(k) + m1.' * (n(k)*W(:,:,k)) * m1 + N/b(k);
    end
    b0(k) = N*S/b0(k);
    
    % ---------------------------------------------------------------------

end


% =========================================================================
% NOT CONSTRAINED
if ~constrained
    
    % ---------------------------------------------------------------------
    %   Gauss-Wishart "precision" parameters
    % ---------------------------------------------------------------------

    for k=1:K
        
        % ---
        % Set up some constants
        sumLogDet = 0;
        sumPsi    = 0;
        Wn        = 0;
        for s=1:S
            [~,~,W,n] = get_po(obj,s);
            sumLogDet = sumLogDet + spm_matcomp('LogDet', W(:,:,k));
            sumPsi    = sumPsi    + spm_prob('DiGamma', n(k)/2, N);
            Wn        = Wn        + n(k)*W(:,:,k);
        end
        sumLogDet = sumLogDet/S;
        sumPsi    = sumPsi/S;
        Wn        = Wn/S;
            
        % -----------------------------------------------------------------
        % Update n0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000
            
            % -------------------------------------------------------------
            % Update W0 (mode, closed-form)
            W0(:,:,k)   = Wn/n0(k);
            LogDetW0(k) = spm_matcomp('LogDet', W0(:,:,k));
            % -------------------------------------------------------------
            
            % ---
            % Objective function
            Eprev = E;
            E = 0.5*S*n0(k)*( LogDetW0(k) - sumLogDet - sumPsi ) ...
                + S*spm_prob('LogGamma', n0(k)/2, N);
            
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = 0.5*S*( LogDetW0(k) - sumLogDet - sumPsi ...
                         + spm_prob('DiGamma', n0(k)/2, N) );
            H = S/4*spm_prob('DiGamma', n0(k)/2, N, 1);

            % ---
            % Update
            n0(k) = max(n0(k) - H\g, N-1+2*eps);
            
        end
        % -----------------------------------------------------------------
    
    end
    
    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    pr.b   = b0;
    pr.m   = m0;
    pr.n   = n0;
    pr.W   = W0;
    pr.ldW = LogDetW0;
    pr.lb  = 0;
    
% =========================================================================
% CONSTRAINED
else

    lb = -inf;
    for em=1:100

        % ---
        % Starting estimate
        if p0 == 0
            p0 = 0;
            V0 = 0;
            for k=1:K
                p0 = p0 + S*n0(k);
                for s=1:S
                    [~,~,W,n] = get_po(obj,s);
                    V0 = V0 + spm_matcomp('Inv', n(k)*W(:,:,k));
                end
            end
            p0 = p0/K;
            V0 = V0/K;
        end

        % -----------------------------------------------------------------
        %   Gauss-Wishart "precision" parameters
        % -----------------------------------------------------------------

        for k=1:K

            % ---
            % Set up some constants
            % > compute sum E[logdet W] and sum psi(nu/2)
            logDetW  = 0;
            psiN     = 0;
            Lambda   = 0;
            for s=1:S
                [~,~,W,n] = get_po(obj,s);
                logDetW = logDetW  + spm_matcomp('Logdet', W(:,:,k));
                psiN    = psiN     + spm_prob('DiGamma', n(k)/2, N);
                Lambda  = Lambda   + n(k)*W(:,:,k);
            end
            logDetW  = logDetW/S;
            psiN = psiN/S;


            % -------------------------------------------------------------
            % Update n0 (mode, Gauss-Newton [convex])
            E = inf;
            for gniter=1:1000

                % ---------------------------------------------------------
                % Update {p,V} for W0 (posterior, closed form)
                p(k)       = p0 + S*n0(k);
                V(:,:,k)   = spm_matcomp('Inv', spm_matcomp('Inv', V0) + Lambda);
                % Useful values
                W0(:,:,k)   = spm_matcomp('Inv', spm_prob('W', 'E', V(:,:,k), p(k)));
                LogDetW0(k) = -spm_prob('W', 'Elogdet', V(:,:,k), p(k));
                % ---------------------------------------------------------

                % ---
                % Objective function
                Eprev = E;
                E = S*n0(k)/2 * (LogDetW0(k) - logDetW - psiN) ...
                    + S*spm_prob('LogGamma', n0(k)/2, N);
                if E == Eprev
                    break;
                end

                % ---
                % Gradient & Hessian
                g = S/2*(LogDetW0(k) - logDetW - psiN + spm_prob('DiGamma', n0(k)/2, N));
                H = S/4 * spm_prob('DiGamma', n0(k)/2, N, 1);

                % ---
                % Update
                n0(k) = max(n0(k) - H\g, N-1+2*eps);
            end
            % ------------------------------------------------------------

        end


        % -----------------------------------------------------------------
        %   Inverse-Wishart parameters
        % -----------------------------------------------------------------

        % ---
        % Set up some constants
        % > compute sum Logdet(psi) and sum psi(m/2)
        sumlogV = 0;
        sumPsi  = 0;
        pV      = 0;
        for k=1:K
            sumlogV = sumlogV + spm_matcomp('LogDet', V(:,:,k));
            sumPsi  = sumPsi  + spm_prob('DiGamma', p(k)/2, N);
            pV      = pV      + p(k)*V(:,:,k);
        end
        sumlogV = sumlogV/K;
        sumPsi  = sumPsi/K;
        pV      = pV/K;


        % -----------------------------------------------------------------
        % Update p0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update V0 (closed-form)
            V0 = pV/p0;
            LogDetV0 = spm_matcomp('LogDet', V0);
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = p0*K/2*( N*LogDetV0 - sumlogV - sumPsi ) + K*spm_prob('LogGamma', p0/2, N);
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = K/2*( LogDetV0 - sumlogV - sumPsi + spm_prob('DiGamma', p0/2, N) );
            H = K/4*spm_prob('DiGamma', p0/2, N, 1);

            % ---
            % Update
            p0 = max(p0 - H\g, N-1+2*eps);

        end
        % -----------------------------------------------------------------
    
        
        % ---
        % Objective function
        lb_prev = lb;
        lb  = 0;
        for k=1:K
            lb  = lb - spm_prob('Wishart', 'kl', V(:,:,k), p(k), V0, p0);
        end
        
        d = abs((lb_prev*(1 + 10*eps) - lb)/lb);
        if 0
            fprintf('%2d | lb = %0.0f | d = %0.5f\n',em,lb,d);  
        end
        
        if d<tol
            break;
        end
        % ---
        
    end % < "EM" loop

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    pr.b   = b0;
    pr.m   = m0;
    pr.n   = n0;
    pr.W   = W0;
    pr.ldW = LogDetW0;
    pr.V   = V;
    pr.p   = p;
    pr.V0  = V0;
    pr.p0  = p0;
    pr.lb  = 0;
    for k=1:K
        pr.lb  = pr.lb - spm_prob('Wishart', 'kl', V(:,:,k), p(k), V0, p0);
    end
    
end

pr.part = part;

% show_all_gaussians(obj, pr)
%==========================================================================


%==========================================================================
function [m,b,W,n] = get_po(obj,s)
m = obj{s}.gmm.po.m;
b = obj{s}.gmm.po.b;
W = obj{s}.gmm.po.W;
n = obj{s}.gmm.po.n;
%==========================================================================


%==========================================================================
function show_all_gaussians(obj, pr, figid)

    if nargin < 3
        figid = 123;
    end
    if figid > 0
        fgauss = figure(figid);
    else
        return
    end
    
    K = numel(pr.b);
    S = numel(obj);
    M = size(pr.W, 1);
    
    if isa(fgauss, 'matlab.ui.Figure')
        
        i = 0;
        for m=1:M
            
            if M > 1
                % ND case
                if m == 1
                    continue;
                end
                ind = [1 m];
            else
                % 1D case
                ind = 1;
            end
            
            
            for k=1:K
                i = i+1;
                subplot(M-1,K,i)
                for s=1:S
                    [mu,~,W,n] = get_po(obj,s);
                    plot_gaussian(mu(ind,k), n(k)*W(ind,ind,k), 'blue', s > 1);
                end
                if isfield(pr, 'V0')
                    h = plot_gaussian(pr.m(ind,k), pr.n(k)*spm_matcomp('Inv', pr.p0*pr.V0), 'black', true);
                    h.LineWidth = 2;
                end
                h = plot_gaussian(pr.m(ind,k), pr.n(k)*pr.W(ind,ind,k), 'red', true);
                h.LineWidth = 2;
            end
        end
        drawnow
    end
%==========================================================================


%==========================================================================
function h = plot_gaussian(mu, Lambda, colour, holdax)

    if nargin < 4
        holdax = false;
        if nargin < 3
            colour = 'b';
        end
    end
    if holdax
        hold on
    end
    
    if numel(mu) > 1
    
        % 2D plot
        
        mu = mu(1:2);
        Sigma2 = spm_matcomp('Inv', Lambda(1:2,1:2));
        Sigma = sqrt(Sigma2);
        [x1,x2] = meshgrid(linspace(mu(1)-3*Sigma(1,1),mu(1)+3*Sigma(1,1),25)', ...
                           linspace(mu(2)-3*Sigma(2,2),mu(2)+3*Sigma(2,2),25)');
        y = mvnpdf([x1(:) x2(:)],mu',Sigma2);
        y = reshape(y, [length(x2) length(x1)]);
        [~,h] = contour(x1,x2,y,1, 'color', colour);
        
    else
        
        % 1D plot
        
        Sigma2 = 1/Lambda;
        Sigma  = sqrt(Sigma2);
        x = linspace(mu-3*Sigma,mu+3*Sigma,25)';
        y = normpdf(x,mu,Sigma2);
        [~,h] = plot(x, y, 'color', colour);
        
    end
    
    if holdax
        hold off
    end
%==========================================================================