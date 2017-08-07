function run_2d(pars)

%==========================================================================
% Initialise algorithm 
%==========================================================================

N = pars.N;
K = pars.K;

samp = pars.samp;
vx   = samp*ones(1,3);

ord = [pars.ord pars.ord pars.ord 0 0 0];

imdir   = pars.imdir;
tempdir = pars.tempdir;

pthmu = pars.pthmu;

debugmode = pars.debugmode;
plotL     = pars.plotL;
showwarp  = pars.showwarp;

docp    = pars.docp;
dobias  = pars.dobias;
doaff   = pars.doaff;
doprior = pars.doprior;
dotpm   = pars.dotpm;

nitmain    = pars.nitmain;
nitcb      = pars.nitcb;
stpaff     = pars.stpaff;
nitaff     = pars.nitaff;
startvelit = pars.startvelit;
tol        = pars.tol;

int_args = pars.int_args;
rparam   = pars.rparam;

% Reg params for affine part (variance)
Raff(1:3,1:3)     = 1e-1*eye(3); % translation
Raff(4:6,4:6)     = 1e-1*eye(3); % rotation
Raff(7:9,7:9)     = 1e3*eye(3);  % scaling
Raff(10:12,10:12) = 1e3*eye(3);  % skew          
 
% Lower bound
L = cell(N,1);
for n=1:N
    L{n} = [-Inf 1];
end
Lglobal = -Inf;

[fig,rndn] = init_fig(N,debugmode,plotL);

tiny = 1e-4;

%--------------------------------------------------------------------------
% Create some directories
%--------------------------------------------------------------------------

if (exist(tempdir,'dir') == 0)
    mkdir(tempdir);
end

tag = ['K-' num2str(K) '_N-' num2str(N) '_vx-' num2str(vx(1)) '_ord-' num2str(ord(1))];

if (exist('results','dir') == 0)
    mkdir('results');
end

f = fullfile('results',tag);
if ~(exist(f,'dir') == 0);
    rmdir(f,'s');
end
mkdir(f);

%--------------------------------------------------------------------------
% Read and process image data
%--------------------------------------------------------------------------

pth = dir(fullfile(tempdir,'preproc','*.nii'));    

if N>numel(pth)
   error('N>numel(pth)');
end

for n=1:N
    pthx{n} = fullfile(tempdir,'preproc',pth(n).name);
end

%--------------------------------------------------------------------------
C = 1;

V = get_V(imdir,tempdir,N);

d     = V{1}.dim;
matim = V{1}.mat;

fmsk = fullfile(tempdir,'msk');
if (exist(fmsk,'dir') == 0)
    mkdir(fmsk);
end

pthmsk = cell(N,1);
for n=1:N 
    Nii  = nifti(pthx{n});
        
    x = Nii.dat(:,:,:);   
    
    % Create mask----------------------------------------------------------
    msk = x~=0 & x~=max(x(:)) & x~=min(x(:)) & isfinite(x);
    
    [~,nam,ext] = fileparts(Nii.dat.fname);
        
    pthmsk{n} = fullfile(fmsk,['msk' nam ext]);
    create_nii(pthmsk{n},msk,matim,'uint8-le','Mask');
end

%--------------------------------------------------------------------------
% Initialise atlas
%--------------------------------------------------------------------------

[mu,matmu] = init_mu(pthx,pthmu,matim,d,K,samp);
lnmu       = log(mu);       

%--------------------------------------------------------------------------
% Initialise cluster parameters
%--------------------------------------------------------------------------

cp = init_cp(pthx,pthmsk,pthmu,mu,K,C); 

%--------------------------------------------------------------------------
% Initialise bias field parameters
%--------------------------------------------------------------------------

chan = init_bf(N,d,vx,C,V,samp);    

%--------------------------------------------------------------------------
% Initialise registration parameters
%--------------------------------------------------------------------------

[pthv,sched,B,affpar] = init_reg(N,d,nitmain,tempdir,matim);

%--------------------------------------------------------------------------
% Set below to true in order to debug registration part
%--------------------------------------------------------------------------

if false   
    docp    = 0;
    dobias  = 0;

    rfolder = '/home/mbrud/Dropbox/PhD/dev/MATLAB/JAs_gems/diffeo_registration/';
    Nmu = nifti(fullfile(rfolder,'Template_4.nii'));
    Nf  = nifti(char(fullfile(rfolder,'rc1mprage0040063.nii'),fullfile(rfolder,'rc2mprage0040063.nii')));

    mu = single(Nmu.dat(:,:,:,:));
    
    lnmu = log(mu);
    
    r  = single(cat(4,Nf(1).dat(:,:,:),Nf(2).dat(:,:,:)));
    r  = cat(4,r,max(1-sum(r,4),0));        
    
    matmu = Nmu.mat;
    matim = Nf.mat;
    
    d = size(Nmu.dat(:,:,:));
    K = size(Nmu.dat(:,:,:,:),4);
    
    M  = matmu(1:3,1:3);
    vx = sqrt(sum(M.^2));
    
    pthmu = 'dummy';

    msk    = ~isfinite(sum(r,4));
    msk    = repmat(msk,1,1,1,K);
    r(msk) = 1/K;
    
    cp(1).w = ones(1,K)/K;
    
    figure(777); clf
end

%==========================================================================
% Start algorithm 
%==========================================================================
 
fprintf('===================================\n')
fprintf('Running algorithm\n')
for itmain=1:nitmain        
            
    muden = 0;
    munum = 0;      
    muv   = 0;
    
    if itmain<startvelit
        dovel = 0;
    elseif itmain==startvelit
        dovel  = 1;
    end
    
    if itmain==stpaff
       nitaff = 1; 
    end
    
    % Decrease regularisation over iterations------------------------------     
    if dovel        
        prm    = [vx rparam*sched(itmain - startvelit + 1)*prod(vx)];
        Greens = spm_shoot_greens('kernel',d,prm);
    else
        prm    = [];
        Greens = [];
    end             
    
    % Mean correcting------------------------------------------------------
    if N>1
        % Affine transform parameters--------------------------------------
        muaffpar = mean(cat(3,affpar{:}),3);
        
        if dovel
            % Velocity fields--------------------------------------------------            
            for n=1:N
               [~,v,~] = load_from_nii(n,pthx,pthv,pthmsk);                         
               muv     = muv + v;
            end
            muv = muv/N;

            clear v
        end
    end
    
    for n=1:N % For each subject
        fprintf('it=%d, n=%d\n',itmain,n); 
        
        r = 0;
        
        % Load data--------------------------------------------------------
        [x,v,msk] = load_from_nii(n,pthx,pthv,pthmsk);                        
%         [x,v] = load_from_nii(n,pthx,pthv,pthmsk);   
%         msk   = ones([prod(d) 1],'logical');
        
        % Warp mu----------------------------------------------------------        
        [phi,~,theta,Jtheta] = make_deformation(v,prm,int_args,Greens,dovel);
        E                    = spm_dexpm(affpar{n},B);
        Affine               = matmu\E*matim;     
        y1                   = affine_transf(Affine,phi);                   
%         wlnmu                = warp(lnmu,y1,ord);       
        wlnmu = sample_mu(lnmu,y1);
                            
        % Get bias field---------------------------------------------------
        bf = get_bf(chan{n},C,d);
        
        % Compute priors of objective--------------------------------------
        [Ebf,Eaff,Evel,Ell] = calc_E(dobias,doaff,dovel,chan{n},itmain,N,pthmu,affpar{n},Raff,v,prm);
      
        
        %------------------------------------------------------------------
        % Cluster and bias part
        %------------------------------------------------------------------
            
        if docp || dobias                                                  
            
            for itcb=1:nitcb
                
                %----------------------------------------------------------
                % Update cluster parameters
                %----------------------------------------------------------   
                    
                if docp                                                    
                    [Ell,r,cp(n)] = update_cp(K,x,bf,wlnmu,cp(n),msk,d);
                    L{n}          = [L{n}; Ell + Ebf + Eaff + Evel 2];
                end

                
                %----------------------------------------------------------
                % Update bias field parameters
                %----------------------------------------------------------
                
                if dobias
                    % Precisions-----------------------------------------------
                    pr = zeros(C,C,K,'single');
                    for k=1:K
                        pr(:,:,k) = inv(cp(n).po.W(:,:,k))/(cp(n).po.nu(k) - C - 1);
                    end

                    for c=1:C
                        if numel(chan{n}(c).T)>0

                            % Get derivatives and Hessian wrt bias parameters--
                            [g,H] = diff_bf(chan{n},bf,r,x,msk,d,K,C,cp(n),pr,c);

                            % Include prior------------------------------------
                            g = g + chan{n}(1).C*chan{n}(c).T(:);
                            H = H + chan{n}(1).C;

                            % Gauss-Newton update------------------------------
                            Update = reshape(H\g,size(chan{n}(c).T));

                            % Do linesearch--------------------------------------------                        
                            oT    = chan{n}(c).T;
                            oEbf  = Ebf;
                            oEll  = Ell;
                            scale = 1;
                            for line_search=1:12
                                chan{n}(c).T = chan{n}(c).T - scale*Update;
                                
                                bf  = get_bf(chan{n},C,d);                                
                                Ebf = -0.5*chan{n}(1).T(:)'*chan{n}(1).C*chan{n}(1).T(:);                                
                                Ell = update_cp(K,x,bf,wlnmu,cp(n),msk,d);

                                if Ell + Ebf + Eaff + Evel > L{n}(end,1)
                                    break;
                                else                                
                                    chan{n}(c).T = oT;
                                    scale    = 0.5*scale;

                                    if line_search==12
                                        % Revert back to old parameters--------
                                        Ebf = oEbf;
                                        Ell = oEll;
                                        bf  = get_bf(chan{n},C,d);
                                    end
                                end 
                            end
                        end
                    end                              

                    L{n} = [L{n}; Ell + Ebf + Eaff + Evel 3];  
                end

                if itcb==nitcb && dobias && docp
                    % Update responsibilities and cluster parameters-------                
                    [Ell,r,cp(n)] = update_cp(K,x,bf,wlnmu,cp(n),msk,d);
                    L{n}          = [L{n};  Ell + Ebf + Eaff + Evel 2];  
                end   
            end    
        end
                
        %------------------------------------------------------------------
        % Registration part
        %------------------------------------------------------------------
        
        if (doaff || dovel) && (~isempty(pthmu) || (itmain>1 && N>1))                                            
                
            %--------------------------------------------------------------
            % Update affine parameters
            %--------------------------------------------------------------                                
            
            if doaff                                                         
                for itaff=1:nitaff
                    
                    r = reshape(r,[d K]);
                
                    if showwarp
                        z = floor((d(3) + 1)/2);
                        figure(777)
                        for k=1:K
                            subplot(2,K,k  ); imagesc(wlnmu(:,:,z,k)'); axis image xy off; colormap(gray);                
                            subplot(2,K,k+K); imagesc(r(:,:,z,k)');    axis image xy off; colormap(gray);
                        end
                        drawnow
                    end
                                            
                    if N>1
                        % Mean correct-------------------------------------
                        affpar{n} = affpar{n} - muaffpar;
                    end
                    
                    % Get derivatives and Hessian wrt affine parameters----
                    [g,H] = diff_aff(affpar{n},r,lnmu,log(cp(n).w),B,matmu,matim,phi,ord,msk);

                    % Include prior----------------------------------------
%                     g = g + Raff*affpar{n};               
%                     H = H + Raff;

                    % Gauss-Newton update----------------------------------
                    Update = H\g; 

                    % Do linesearch----------------------------------------    
%                     Pmu    = weighted_softmax(wlnmu,log(cp.w),msk);
%                     Ell    = sum(sum(rmsk(r,msk,d,K).*log(Pmu),2),1);
                    
                    oEaff   = Eaff;
                    oEll    = Ell;
                    oaffpar = affpar{n};
                    scale   = 1;
                    for line_search=1:12
                        affpar{n} = affpar{n} - scale*Update;
                        
                        E      = spm_dexpm(affpar{n},B);
                        Affine = matmu\E*matim; 
                        y1     = affine_transf(Affine,phi);                    
%                         wlnmu  = warp(lnmu,y1,ord);                     
                        wlnmu = sample_mu(lnmu,y1);
        
%                         Pmu   = weighted_softmax(wlnmu,log(cp.w),msk);
%                         Ell   = sum(sum(rmsk(r,msk,d,K).*log(Pmu),2),1)
                        
                        Ell  = update_cp(K,x,bf,wlnmu,cp(n),msk,d); 
                        Eaff = 0;%-0.5*affpar{n}'*Raff*affpar{n}; 

                        if Ell + Ebf + Eaff + Evel > L{n}(end,1)
                            if showwarp
                                fprintf('Affine converged (n=%d,it=%d)\n',n,line_search); 

                                z = floor((d(3) + 1)/2);
                                figure(777)
                                for k=1:K
                                    subplot(2,K,k  ); imagesc(wlnmu(:,:,z,k)'); axis image xy off; colormap(gray);
                                    subplot(2,K,k+K); imagesc(r(:,:,z,k)');     axis image xy off; colormap(gray);
                                end
                                drawnow
                            end                        

                            [Ell,r,cp(n)] = update_cp(K,x,bf,wlnmu,cp(n),msk,d);       
                            L{n}          = [L{n}; Ell + Ebf + Eaff + Evel 4];  
                        
                            break;
                        else
                            scale     = 0.5*scale;
                            affpar{n} = oaffpar;

                            if line_search==12
                                % Revert back to old parameters------------
                                E      = spm_dexpm(affpar{n},B);
                                Affine = matmu\E*matim; 
                                y1     = affine_transf(Affine,phi);
%                                 wlnmu  = warp(lnmu,y1,ord); 
                                wlnmu = sample_mu(lnmu,y1);
        
                                Eaff = oEaff;
                                Ell  = oEll;
                                
                                L{n} = [L{n}; Ell + Ebf + Eaff + Evel 4];
                            end
                        end
                    end
                end                                
            end
            
            %--------------------------------------------------------------
            % Update velocity parameters
            %--------------------------------------------------------------
            
            if dovel                                
                
                r = reshape(r,[d K]);
                
                if showwarp
                    z = floor((d(3) + 1)/2);
                    figure(777)
                    for k=1:K
                        subplot(2,K,k  ); imagesc(wlnmu(:,:,z,k)'); axis image xy off; colormap(gray);                
                        subplot(2,K,k+K); imagesc(r(:,:,z,k)');    axis image xy off; colormap(gray);
                    end
                    drawnow
                end
                
                if N>1
                    % Mean correct-----------------------------------------
                    v = v - muv;
                end
                    
                % Get derivatives and Hessian wrt velocity parameters------
                [g,H] = diff_vel(r,lnmu,log(cp(n).w),Affine,phi,int_args,ord,msk);

                % Include prior--------------------------------------------
                u = spm_diffeo('vel2mom',v,prm);
                g = g + u;

                % Gauss-Newton update--------------------------------------
                Update = spm_diffeo('fmg',H,g,[prm 3 2]);

                % Do linesearch--------------------------------------------
%                 Pmu   = weighted_softmax(wlnmu,log(cp.w),msk);
%                 Ell   = sum(sum(rmsk(r,msk,d,K).*log(Pmu)));
                
                oEll  = Ell;
                oEvel = Evel;
                ov    = v;
                scale = 1;
                for line_search=1:12
                    v = v - scale*Update;
                    
                    [phi,~,theta,Jtheta] = make_deformation(v,prm,int_args,Greens,dovel);
                    y1                   = affine_transf(Affine,phi);                    
%                     wlnmu                = warp(lnmu,y1,ord); 
                    wlnmu = sample_mu(lnmu,y1);
                              
%                     Pmu = weighted_softmax(wlnmu,log(cp.w),msk);
%                     Ell = sum(sum(rmsk(r,msk,d,K).*log(Pmu)));

                    Ell  = update_cp(K,x,bf,wlnmu,cp(n),msk,d);
                    u    = spm_diffeo('vel2mom',v,prm);
                    Evel = -0.5*sum(sum(sum(sum(u.*v))));
                    
                    if Ell + Ebf + Eaff + Evel  > L{n}(end,1) 
                        if showwarp
                            fprintf('Nonlinear converged (n=%d,it=%d)\n',n,line_search); 
                            
                            z = floor((d(3) + 1)/2);
                            figure(777)
                            for k=1:K
                                subplot(2,K,k  ); imagesc(wlnmu(:,:,z,k)'); axis image xy off; colormap(gray);
                                subplot(2,K,k+K); imagesc(r(:,:,z,k)');     axis image xy off; colormap(gray);
                            end
                            drawnow
                        end
                        
                        [Ell,r,cp(n)] = update_cp(K,x,bf,wlnmu,cp(n),msk,d);       
                        L{n}          = [L{n}; Ell + Ebf + Eaff + Evel 5];  
                            
                        break;
                    else
                        scale = 0.5*scale;
                        v     = ov;
                        
                        if line_search==12
                            % Revert back to old parameters----------------
                            [~,~,theta,Jtheta] = make_deformation(v,prm,int_args,Greens,dovel);
                            
                            Evel = oEvel;                            
                            Ell  = oEll;
                            
                            L{n} = [L{n}; Ell + Ebf + Eaff + Evel 5];  
                        end
                    end
                end
            end     

        end    
        
        %------------------------------------------------------------------
        % Update statistics for TPM
        %------------------------------------------------------------------
        
        if dotpm && N>1                                
            % Compute multinomial part
            wmu  = bsxfun(@times,reshape(mu,[prod(d),K]),cp(n).w);
            mnom = bsxfun(@rdivide,wmu,sum(wmu,2));
            
            % Compute likelihood part
            [~,~,~,lik] = update_cp(K,x,bf,wlnmu,cp(n),msk,d); 
            
            % Warp likelihood to template space
            E      = spm_dexpm(affpar{n},B);
            Affine = matmu\E*matim; 
            y1     = spm_diffeo('comp',theta,affine_transf(inv(Affine),identity(d)));
            
            y1(:,:,:,3) = 1;
            
            lik    = warp(reshape(lik,[d K]),y1,ord);   
            
            % Multiply warped likelihood w. determinants of inv(Affine) and Jtheta
            lik = bsxfun(@times,lik,spm_diffeo('det',Jtheta)*det(inv(Affine)));                                            
            lik = reshape(lik,[prod(d) K]);
            
            % Compute numerator part of mu update for subject n
            tmp             = mnom.*lik;    
            tmp             = bsxfun(@rdivide,tmp,sum(tmp,2));
            tmp(isnan(tmp)) = 0;                        
            munum           = munum + tmp; % Sum over subjects
        
            % Compute denomenator part of mu update for subject n
            tmp             = 1./sum(wmu,2);    
            tmp             = bsxfun(@times,cp(n).w,repmat(tmp,1,K));
            tmp(isnan(tmp)) = 0;                              
            muden           = muden + tmp; % Sum over subjects
        end
        
        % Save data--------------------------------------------------------
        save_to_nii(n,pthv,v);
        
        fprintf('cp(%i).w: %s\n',n,sprintf('%3.3f ',cp(n).w));
    end
    
    %----------------------------------------------------------------------
    % Normalize TPMs
    %----------------------------------------------------------------------            
    
    if dotpm && N>1
%         % "better way" (John hopes)----------------------------------------
%         mu = munum./muden; 
%         mu = bsxfun(@rdivide,mu,sum(mu,2)); 
        
        % "naive way"------------------------------------------------------
        mu = bsxfun(@rdivide,munum,sum(munum,2)); 
        
        % Make sure that mu is a valid probability distribution------------
        msk     = ~isfinite(sum(mu,2)) | sum(mu,2) == 0;
        msk     = repmat(msk,1,K);
        mu(msk) = 1/K;                
        
        if sum(~isfinite(mu(:))), 
            warning('sum(~isfinite(mu(:)))'); 
        end        
                
        mu   = reshape(mu,[d K]);
        lnmu = log(mu + tiny);
    
        if sum(~isfinite(lnmu(:))), 
            warning('sum(~isfinite(lnmu(:)))'); 
        end        
        
        % Compute objective------------------------------------------------                
        for n=1:N
            [x,v,msk] = load_from_nii(n,pthx,pthv,pthmsk);
            
            [Ebf,Eaff,Evel] = calc_E(dobias,doaff,dovel,chan{n},itmain,N,pthmu,affpar{n},Raff,v,prm);  
            
            phi = make_deformation(v,prm,int_args,Greens,dovel);
            
            E      = spm_dexpm(affpar{n},B);
            Affine = matmu\E*matim;                           
            y1     = affine_transf(Affine,phi);                    
%             wlnmu  = warp(lnmu,y1,ord); 
            wlnmu = sample_mu(lnmu,y1);
                          
            bf = get_bf(chan{n},C,d);
            
            Ell = update_cp(K,x,bf,wlnmu,cp(n),msk,d);  
            
            L{n} = [L{n}; Ell + Ebf + Eaff + Evel 6];
                                
            % Update responsibilities and cluster parameters---------------
            [Ell,~,cp(n)] = update_cp(K,x,bf,wlnmu,cp(n),msk,d);       
            L{n}          = [L{n}; Ell + Ebf + Eaff + Evel 2];              
        end
    end
        
    %----------------------------------------------------------------------
    % Update intensity priors
    %----------------------------------------------------------------------    
    
    if doprior && N>1 
        oLpr = cat(3,L{:});
        oLpr = sum(oLpr(end,1,:));
        ocp  = cp;
        
        pr = update_pr(cp);
        for n=1:N
            cp(n).pr = pr;
        end
        
        % Compute objective------------------------------------------------  
        L1 = zeros(N,2);
        L2 = zeros(N,2);        
        for n=1:N
            [x,v,msk] = load_from_nii(n,pthx,pthv,pthmsk);
            
            [Ebf,Eaff,Evel] = calc_E(dobias,doaff,dovel,chan{n},itmain,N,pthmu,affpar{n},Raff,v,prm);     
            
            phi = make_deformation(v,prm,int_args,Greens,dovel);
            
            E      = spm_dexpm(affpar{n},B);
            Affine = matmu\E*matim;                           
            y1     = affine_transf(Affine,phi);                    
%             wlnmu  = warp(lnmu,y1,ord); 
            wlnmu = sample_mu(lnmu,y1);
                        
            bf = get_bf(chan{n},C,d);
            
            Ell = update_cp(K,x,bf,wlnmu,cp(n),msk,d);  
            
            L1(n,:) = [Ell + Ebf + Eaff + Evel 7];
            
            % Update responsibilities and cluster parameters---------------
            [Ell,~,cp(n)] = update_cp(K,x,bf,wlnmu,cp(n),msk,d);       
            L2(n,:)       = [Ell + Ebf + Eaff + Evel 2];              
        end    
        
        nLpr = sum(L1(:,1));
        
        if nLpr<oLpr
            cp = ocp;
        else
            for n=1:N
                L{n} = [L{n}; L1(n,:)];
                L{n} = [L{n}; L2(n,:)];
            end
        end
    end    
                  
    %----------------------------------------------------------------------
    % Plot, check stopping criteria, save, etc.
    %----------------------------------------------------------------------
    
    % Display results------------------------------------------------------    
    if debugmode
        visualize_progress(fig,d,K,rndn,mu,pthv,C,chan,affpar,prm,int_args,Greens,dovel,matmu,matim,ord,B,pthx,pthmsk,cp);
    end
        
    % Save most recent TPM-------------------------------------------------
    write_mu(tag,mu,matmu);
        
    
    % Plot lower bound-----------------------------------------------------
    if plotL
        plot_L(L,itmain,fig);  
    end
    
    % Check convergence----------------------------------------------------
    Lglobal = get_Lglobal(Lglobal,L);                  
    
    fprintf('%d %g %g\n',itmain,Lglobal(end),abs((Lglobal(end) - Lglobal(end - 1))/Lglobal(end)));
    
    if abs((Lglobal(end) - Lglobal(end - 1))/Lglobal(end))<tol && itmain>20
        fprintf('Algorithm converged in %d iterations.\n',itmain)                        
        break;
    end
end
%==========================================================================

%==========================================================================
function write_mu(tag,mu,matmu)
f     = fullfile('results',tag);
fname = fullfile(f,[tag '__tpm.nii']);
create_nii(fname,mu,matmu,'float32','TPM');   
%==========================================================================

%=======================================================================        
function [Ebf,Eaff,Evel,Ell] = calc_E(dobias,doaff,dovel,chan,itmain,N,pthmu,affpar,Raff,v,prm)
Ebf  = 0;
Eaff = 0;
Evel = 0;
Ell  = 0;

if dobias                                  
    Ebf = -0.5*chan(1).T(:)'*chan(1).C*chan(1).T(:);
end

if doaff && (~isempty(pthmu) || (itmain>1 && N>1))   
    Eaff = 0;%-0.5*affpar'*Raff*affpar;   
end
if dovel && (~isempty(pthmu) || (itmain>1 && N>1))   
    Evel = -0.5*sum(sum(sum(sum(spm_diffeo('vel2mom',v,prm).*v))));
end     
%=======================================================================

%==========================================================================
function [x,v,msk] = load_from_nii(n,pthx,pthv,pthmsk)
x = nifti(pthx{n});
x = single(x.dat(:,:,:));

v = nifti(pthv{n});
v = single(v.dat(:,:,:,:));

msk = nifti(pthmsk{n});
msk = logical(reshape(msk.dat(:,:,:),[],1));
%==========================================================================

%==========================================================================
function save_to_nii(n,pthv,v)
Nii              = nifti(pthv{n});
Nii.dat(:,:,:,:) = v;
%==========================================================================

%==========================================================================
function visualize_progress(fig,d,K,rndn,mu,pthv,C,chan,affpar,prm,int_args,Greens,dovel,matmu,matim,ord,B,pthx,pthmsk,cp)
z = floor((d(3)+1)/2);   
n = rndn(end);

% Velocity fields
set(0,'CurrentFigure',fig{2});

cnt = 1;
for f1=1:numel(rndn)
    Nii = nifti(pthv(rndn(f1)));
    v   = Nii.dat(:,:,:,:);

    for f2=1:3
        subplot(numel(rndn),3,cnt)
        imagesc(v(:,:,z,f2)'); axis image xy off; title(['v' num2str(f2) ', n=' num2str(rndn(f1))]); colormap(gray); colorbar;
        cnt = cnt + 1;
    end
end
drawnow

% Bias fields
set(0,'CurrentFigure',fig{3});

cnt = 1;
for f1=1:numel(rndn)
    bf = get_bf(chan{rndn(f1)},C,d);
    bf = reshape(bf,[d C]);

    bf1 = bf(:,:,floor((size(bf,3) + 1)/2));
    bf2 = squeeze(bf(floor((size(bf,1) + 1)/2),:,:));
    bf3 = squeeze(bf(:,floor((size(bf,2) + 1)/2),:));  

    subplot(numel(rndn),3,cnt)
    imagesc(bf1'); axis image xy off; title(['bf' num2str(1) ', n=' num2str(rndn(f1))]); colormap(parula); colorbar;
    cnt = cnt + 1;

    subplot(numel(rndn),3,cnt)
    imagesc(bf2'); axis image xy off; title(['bf' num2str(2) ', n=' num2str(rndn(f1))]); colormap(parula); colorbar;
    cnt = cnt + 1;

    subplot(numel(rndn),3,cnt)
    imagesc(bf3'); axis image xy off; title(['bf' num2str(3) ', n=' num2str(rndn(f1))]); colormap(parula); colorbar;
    cnt = cnt + 1;    
end    
drawnow

% Warped atlas
set(0,'CurrentFigure',fig{4});

[x,v,msk] = load_from_nii(n,pthx,pthv,pthmsk);

phi    = make_deformation(v,prm,int_args,Greens,dovel);
E      = spm_dexpm(affpar{n},B);
Affine = matmu\E*matim;
y1     = affine_transf(Affine,phi);                    
% wmu    = warp(mu,y1,ord);
wmu = sample_mu(mu,y1);
        
for k=1:K
    subplot(2,K,k);
    imagesc(mu(:,:,z,k)'); axis image xy off; title(['mu, k=' num2str(k)]); colormap(gray);            
end
for k=(K + 1):2*K
    subplot(2,K,k);
    imagesc(wmu(:,:,z,k - K)'); axis image xy off; title(['wmu, n=' num2str(n) ', k=' num2str(k - K)]); colormap(gray);            
end
drawnow

% Resps
set(0,'CurrentFigure',fig{5});    

% wlnmu = warp(log(mu),y1,ord);
% wlnmu = sample_mu(log(mu),y1);
[~,r] = update_cp(K,x,bf,log(mu + 1e-4),cp(n),msk,d);
r     = reshape(r,[d K]);

F  = K;
F1 = floor(sqrt(F));
F2 = ceil(F/F1);      
for f=1:F
    subplot(F1,F2,f);
    imagesc(r(:,:,z,f)'); axis image xy off; title(['r, n=' num2str(n) ', k=' num2str(f)]); colormap(gray);
end  
drawnow
%==========================================================================

%==========================================================================
function Lglobal = get_Lglobal(Lglobal,L)
sL = zeros(size(L{1},1),1);
for n=1:numel(L)
    sL = sL + L{n}(:,1);
end
    
Lglobal(end + 1) = sL(end);

if Lglobal(end) - Lglobal(end - 1) < 0
    fprintf(2,'Global LB decreased by %1.2f\n',Lglobal(end) - Lglobal(end - 1));
end
%==========================================================================

%==========================================================================
function plot_L(L,iter,fig,ix)
if nargin<4, ix = 1; end

sL = zeros(size(L{1},1),1);
for n=1:numel(L)
    sL = sL + L{n}(:,1);
end

n      = 1;
colmap = {'m.','g.','y.','b.','r.','c.','k.'};

set(0,'CurrentFigure',fig{1});

sL  = sL(2:end);
slb = L{n}(2:end,1);
col = L{n}(2:end,2);

subplot(121)

hold on
plot(1:numel(sL),sL,'k-','LineWidth',1);
for i=1:numel(sL)
    plot(i,sL(i,ix),colmap{col(i)},'markersize',10);
end
hold off
xlabel('Convergence (global)')    
title(['iter=' num2str(iter)])

subplot(122)

hold on
plot(1:numel(slb),slb,'k-','LineWidth',1);
for i=1:numel(slb)
    plot(i,slb(i,ix),colmap{col(i)},'markersize',10);
end
hold off
xlabel(['Convergence (n=' num2str(n) ')'])    
title(['iter=' num2str(iter)])

drawnow          
%==========================================================================

%==========================================================================
function [fig,rndn] = init_fig(N,debugmode,plotL)
fig = cell(5,1);

if debugmode && plotL
    for i=1:5
        fig{i} = figure(666 + i); 
        clf(fig{i});     
    end    
elseif debugmode && ~plotL
    for i=1:5
        if i==1
            close(figure(666 + i));
        else
            fig{i} = figure(666 + i); 
            clf(fig{i});
        end
    end    
elseif ~debugmode && plotL
    for i=1:5
        if i==1
            fig{i} = figure(666 + i); 
            clf(fig{i});
        else
           close(figure(666 + i));
        end
    end
else
    for i=1:5
        close(figure(666 + i));
    end
end
drawnow;

rndn = datasample(1:N,min(N,4),'Replace',false);
%==========================================================================

%=======================================================================
function V = get_V(imdir,tempdir,N)
pth = dir(fullfile(imdir,'*.nii'));  
if N>numel(pth)
   error('N>numel(pth)');
end
V = cell(N,2);
for n=1:N
    V{n,1} = spm_vol(fullfile(imdir,pth(n).name));
end
for n=1:N
    V{n,2} = spm_vol(fullfile(tempdir,'preproc',['im' num2str(n) '.nii']));
end
%=======================================================================