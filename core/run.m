function run(pars)

%==========================================================================
% Initialise algorithm 
%==========================================================================

N = pars.N;
K = pars.K;
C = pars.C;

ct = pars.ct;

if pars.runpar
    runpar = get_nbr_of_cores; 
else
    runpar = 0;
end

samp = pars.samp;
vx   = samp*ones(1,3);

bs = pars.bs;

figix      = pars.figix;
debuglevel = pars.debuglevel;
if runpar && debuglevel==4
    debuglevel = 3;    
end

imdir   = pars.imdir;
tempdir = pars.tempdir;

preproc = pars.preproc;

pthmu = pars.pthmu;

do = pars.do;

nitmain    = pars.nitmain;
nitcpbf    = pars.nitcpbf;
nitcp0     = pars.nitcp0;
nitsrtw    = pars.nitsrtw;
nitstrtreg = pars.nitstrtreg;
tol        = pars.tol;

bflam  = pars.bflam;
bffwhm = pars.bffwhm;

abasis = pars.abasis;
alam   = pars.alam;     

vlam = pars.vlam;

mufwhm = samp*pars.mufwhm;

%--------------------------------------------------------------------------
% Read and preprocess image data
%--------------------------------------------------------------------------

if ~preproc.imload 
    % Read and preprocess images from folder imdir------------------------- 

    addpath(genpath('./preproc'))
    
    [pthf,C,N,V] = get_pth2f(imdir,N,C);        
        
    if (exist(tempdir,'dir') == 0)
        mkdir(tempdir);
    end

    if numel(dir(tempdir)) > 2
        rmdir(fullfile(tempdir,'*'),'s'); % Clear tempdir
    end
    
    % Preprocess images (reset origin, align to MNI, etc) -----------------
    pthf = preproc_im(pthf,preproc,tempdir,runpar);
else
    % Load already preprocessed images-------------------------------------
    
    [pthf,C,N,V] = get_pth2f(fullfile(tempdir,'preproc'),N,C);  
end

%--------------------------------------------------------------------------
% Warp images to same size
%--------------------------------------------------------------------------

[pthf,Mf,d] = warp_eq_size(pthf,samp,tempdir,ct);

%--------------------------------------------------------------------------
% Initialise atlas
%--------------------------------------------------------------------------

[lnmu,Mmu] = init_mu(pthf,pthmu,Mf,d,K,samp);

%--------------------------------------------------------------------------
% Initialise cluster struct
%--------------------------------------------------------------------------

cp = init_cp(pthf,pthmu,lnmu,K,d,runpar); 

%--------------------------------------------------------------------------
% Initialise bias field struct
%--------------------------------------------------------------------------

chan = init_bf(N,C,V,d,bflam,bffwhm);    

%--------------------------------------------------------------------------
% Initialise registration parameters
%--------------------------------------------------------------------------

[pthv,sched,B,a,R,prm0,prm,int_args,vconv,do] = init_reg(N,d,nitmain,tempdir,Mf,abasis,alam,vx,vlam,do);

%--------------------------------------------------------------------------
% Miscellaneous (results directory, lower bound, visualisation, etc)
%--------------------------------------------------------------------------

if do.writemu
    tag = ['K-' num2str(K) '_N-' num2str(N) '_vx-' num2str(vx(1)) '_deg-' num2str(bs(1))];

    if (exist('results','dir') == 0)
        mkdir('results');
    end

    folder = fullfile('results',tag);
    if ~(exist(folder,'dir') == 0)
        rmdir(folder,'s');
    end
    mkdir(folder);
end

[L,gL] = init_L(N);

[fig,rndn,zix] = init_fig(N,debuglevel,d,figix);
fig5           = fig{5};

if true
    % Show image data (for debugging)
    figure(1 + 5*(figix - 1) + 5);
    NN = min(N,32);
    N1 = floor(sqrt(NN));
    N2 = ceil(NN/N1);      
    for n=1:NN
        subplot(N1,N2,n);
        
        Nii = nifti(pthf{n,1});
        img = Nii.dat(:,:,zix);
        imagesc(img'); axis image xy off; colormap(gray); colorbar;
    end  
    drawnow
end

%==========================================================================
% Start algorithm 
%==========================================================================
 
fprintf('===================================\n')
fprintf('Running algorithm\n')
for itmain=1:nitmain                          

    if itmain==2
        % Most of the objfun improvements are in the first iteration,
        % so show only improvements after this, as they are more clearly visible.
        [L,gL] = init_L(N);
    end
    
    % Increase the number of MoG iterations over time, so that the
    % intensity priors burn in and so that there is time for the 
    % template knowledge to propagate between subjects.
    if itmain>numel(nitcp0)
        nitcp = nitcp0(end);
    else
        nitcp = nitcp0(itmain);
    end
    
    if itmain==nitsrtw && do.w
        % Start weight update
        for n=1:N
            cp(n).dow = 1;
        end
    end
    
    if itmain==nitstrtreg(3) && do.v0
        % Large deformation model
        int_args = 8;
        do.v     = 1;
    end
    if itmain==nitstrtreg(2) && do.v0
        % Small deformation model
        int_args = 1;
        do.v     = 1;
    end
    if itmain==nitstrtreg(1) && do.a0
        % Affine
        do.a = 1;
    end
       
    % Depending on which component of the model is updated different Green
    % functions as well as regularisation parameters need to be used. This
    % is so that the objective functions converges properly.
    Greens0 = [];
    Greens  = [];                   
    if int_args>1
        % Decrease velocity field regularisation over iterations
        prm0 = [vx vlam*sched(itmain - nitstrtreg(3) + 1)];                 
        prm  = [vx vlam*sched(itmain - nitstrtreg(3) + 2)]; 

        if itmain>nitstrtreg(3)
            Greens0 = spm_shoot_greens('kernel',d,prm0);    
        end
        Greens = spm_shoot_greens('kernel',d,prm);    
    end   
    
    % Mean correct affine transform parameters
    if N>1 && do.a && do.amu, amu = mean(cat(3,a{:}),3);
    else                      amu = 0; end 
    
    muden = 0; munum = 0; % Used for computing new template 
        
    parfor (n=1:N,runpar) % For each subject
%     for n=1:N
        fprintf('it=%d, n=%d\n',itmain,n);                   
    
        r = 0; Lz = 0; % Needs to be defined to avoid warning message from parfor...
        
        % Load data--------------------------------------------------------
        [f,v,msk] = load_from_nii(n,pthf,pthv,d,C,ct);            
        
        % Warp template to subject-----------------------------------------  
        [phi,~,theta,Jtheta] = make_deformation(v,prm0,int_args,Greens0);
        E                    = spm_dexpm(a{n},B);
        Affine               = Mmu\E*Mf;     
        y1                   = affine_transf(Affine,phi);                   
        wlnmu                = warp(lnmu,y1,bs); y1 = [];        
        
        % Re-generate bias field-------------------------------------------
        bf = get_bf(chan{n},d);
        
        % Objective function-----------------------------------------------
        [~,Lbf,La,Lv] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm0,msk,false);

        %------------------------------------------------------------------
        % Cluster and bias part
        %------------------------------------------------------------------                                                          
            
        for itcpbf=1:nitcpbf

            %----------------------------------------------------------
            % Update cluster parameters
            %----------------------------------------------------------   
                    
            for itcp=1:nitcp
                [Lz,r,cp(n)] = update_cp(f,bf,wlnmu,cp(n),msk); r = reshape(r,[d K]);
                L{n}         = [L{n}; Lz + Lbf + La + Lv, 2];
            end
            
            %----------------------------------------------------------
            % Update bias field parameters (as in spm_preproc8)
            %----------------------------------------------------------

            if do.bf
            
                % Precisions-----------------------------------------------
                cpLam = get_cpLam(cp(n));        

                oL = L{n}(end,1);
                for c=1:C                        
                    f  = reshape(f,[d C]);
                    bf = reshape(bf,[d C]);     
                    
                    d3 = numel(chan{n}(c).T);
                    if d3>0
                        
                        % Gradient and Hessian------------------------------
                        H = zeros(d3,d3);
                        g = zeros(d3,1);  

                        for z=1:d(3)

                            mskz = msk(:,:,z,c);
                            nm   = nnz(mskz);

                            if nm==0, continue; end
                            
                            cr = cell(N,1);
                            for c1=1:C
                                fz  = f(:,:,z,c1);  fz  = fz(mskz);
                                bfz = bf(:,:,z,c1); bfz = bfz(mskz);
                                
                                cr{c1} = double(fz).*double(bfz); 
                            end

                            w1 = zeros(nm,1);
                            w2 = zeros(nm,1);
                            for k=1:K
                                rk  = r(:,:,z,k);
                                rk  = rk(mskz);
                                w0  = zeros(nm,1);
                                for c1=1:C
                                    w0 = w0 + cpLam(c1,c,k)*(cp(n).po.m(c1,k) - cr{c1});
                                end
                                w1  = w1 + rk.*w0;
                                w2  = w2 + rk*cpLam(c,c,k);
                            end
                            wt1       = zeros(d(1:2));
                            wt1(mskz) = -(1 + cr{c}.*w1); % US eq. 34 (gradient)
                            wt2       = zeros(d(1:2));
                            wt2(mskz) = cr{c}.*cr{c}.*w2 + 1; % Simplified Hessian of US eq. 34
                            cr = [];

                            b3 = chan{n}(c).B3(z,:)';
                            g  = g  + kron(b3,spm_krutil(wt1,chan{n}(c).B1,chan{n}(c).B2,0));
                            H  = H + kron(b3*b3',spm_krutil(wt2,chan{n}(c).B1,chan{n}(c).B2,1));
                            wt1 = []; wt2 = []; b3 = [];
                        end

                        bfLam = chan{n}(c).C; % Inverse covariance of priors                            

                        % Add prior----------------------------------------
                        g = g + bfLam*chan{n}(c).T(:);               
                        H = H + bfLam;

                        % Gauss-Newton update------------------------------
                        Update = H\g; g = []; H = [];                            
                        Update = reshape(Update,size(chan{n}(c).T));
                        
                        % Line-search--------------------------------------                        
                        f = reshape(f,[prod(d) C]);
                        
                        scale = 1.0;
                        oLz   = Lz;
                        oLbf  = Lbf;
                        oT    = chan{n}(c).T;
                        for line_search=1:12
                            chan{n}(c).T = chan{n}(c).T - scale*Update; % Backtrack if necessary

                            % Re-generate bias field-----------------------
                            bf = get_bf(chan{n},d);  

                            % Objective function---------------------------
                            [Lz,Lbf] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm0,msk);

                            if Lz + Lbf + La + Lv >= oL
                                oL = Lz + Lbf + La + Lv;
                                break;
                            else                                
                                chan{n}(c).T = oT;
                                scale        = 0.5*scale;

                                if line_search==12
                                    % Revert back to old parameters--------
                                    Lbf = oLbf;
                                    Lz  = oLz;
                                    bf  = get_bf(chan{n},d);
                                end
                            end 
                        end                        
                        oT = []; Update = [];
                    end                                                                
                end                
                cpLam = [];
                
                L{n} = [L{n}; Lz + Lbf + La + Lv, 3]; 
                
                if itcpbf==nitcpbf
                    % Update responsibilities and cluster parameters-------                
                    [Lz,r,cp(n)] = update_cp(f,bf,wlnmu,cp(n),msk); r = reshape(r,[d K]);
                    L{n}         = [L{n};  Lz + Lbf + La + Lv, 2];  
                end 
            end 
        end   
             
        if debuglevel>3
            set(0,'CurrentFigure',fig5);                                              
            show_mu(fig5,K,wlnmu,r,d,zix);
        end
            
        %------------------------------------------------------------------
        % Registration part
        %------------------------------------------------------------------
        
        if ~isempty(pthmu) || (itmain>1 && N>1)                                                            
            
            %--------------------------------------------------------------
            % Update affine parameters
            %--------------------------------------------------------------                                
            
            if do.a            
                
                if do.amu 
                    % Mean correct affine parameters-------------------------------
                    a{n} = a{n} - amu;
                end
                
                % Gradient and Hessian-------------------------------------
                [g,H] = diff_a(a{n},r,lnmu,cp(n).lnw,B,Mmu,Mf,phi);

                % Add prior----------------------------------------------
                g = g + R*a{n};               
                H = H + R;
                    
                % Gauss-Newton update--------------------------------------
                Update = H\g; g = []; H = []; 

                % Line-search----------------------------------------------
                scale = 1.0;
                oLa   = La;
                oLz   = Lz;
                oa    = a{n};
                for line_search=1:12
                    a{n} = a{n} - scale*Update; % Backtrack if necessary

                    % Warp template to subject-----------------------------
                    E      = spm_dexpm(a{n},B);
                    Affine = Mmu\E*Mf; 
                    y1     = affine_transf(Affine,phi);                    
                    wlnmu  = warp(lnmu,y1,bs); y1 = [];                   

                    % Objective function-----------------------------------
                    [Lz,Lbf,La] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm0,msk);

                    if Lz + Lbf + La + Lv >= L{n}(end,1)

                        if debuglevel>3
                            fprintf('Affine converged (n=%d,it=%d)\n',n,line_search); 
                            show_mu(fig5,K,wlnmu,r,d,zix);
                        end                        

                        break;
                    else
                        scale = 0.5*scale;
                        a{n}  = oa;

                        if line_search==12
                            % Revert back to old parameters----------------
                            E      = spm_dexpm(a{n},B);
                            Affine = Mmu\E*Mf; 
                            y1     = affine_transf(Affine,phi);
                            wlnmu  = warp(lnmu,y1,bs); y1 = [];           

                            La = oLa;
                            Lz = oLz; 
                        end
                    end
                end  

                L{n} = [L{n}; Lz + Lbf + La + Lv, 4];
                
                [Lz,r,cp(n)] = update_cp(f,bf,wlnmu,cp(n),msk); r = reshape(r,[d K]);       
                L{n}         = [L{n}; Lz + Lbf + La + Lv, 2];

                Update = [];                              
            end
            
            %--------------------------------------------------------------
            % Update velocity parameters
            %--------------------------------------------------------------
            
            if do.v                                        
                
                % Gradient and Hessian-------------------------------------
                [g,H] = diff_v(r,lnmu,cp(n).lnw,Affine,phi,int_args);

                % Add prior------------------------------------------------
                u = spm_diffeo('vel2mom',v,prm);
                g = g + u; u = [];

                % Gauss-Newton update--------------------------------------
                Update = spm_diffeo('fmg',H,g,[prm 3 2]); g = []; H = [];

                % Line-search----------------------------------------------
                scale = 1.0;
                oLz   = Lz;
                oLv   = Lv;
                ov    = v;
                for line_search=1:12
                    v = v - scale*Update; % Backtrack if necessary
                    if d(3)==1
                        v(:,:,:,3) = 0;
                    end
                    
                    % Warp template to subject-----------------------------
                    [phi,~,theta,Jtheta] = make_deformation(v,prm,int_args,Greens);
                    y1                   = affine_transf(Affine,phi);                    
                    wlnmu                = warp(lnmu,y1,bs); y1 = [];          
                              
                    % Objective function-----------------------------------
                    [Lz,Lbf,La,Lv] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm,msk);
                    
                    if Lz + Lbf + La + Lv >= L{n}(end,1) 
                        
                        if debuglevel>3
                            fprintf('Nonlinear converged (n=%d,it=%d)\n',n,line_search); 
                            show_mu(fig5,K,wlnmu,r,d,zix);                            
                        end                        
                            
                        vconv(n) = true;
                        break;
                    else
                        scale = 0.5*scale;
                        v     = ov;
                        
                        if line_search==12
                            % Revert back to old parameters----------------
                            [phi,~,theta,Jtheta] = make_deformation(v,prm0,int_args,Greens0);
                            
                            E      = spm_dexpm(a{n},B);
                            Affine = Mmu\E*Mf; 
                            y1     = affine_transf(Affine,phi);
                            wlnmu  = warp(lnmu,y1,bs); y1 = [];    
                                
                            Lv = oLv;                            
                            Lz = oLz;      
                            
                            vconv(n) = false;
                        end
                    end
                end
                
                L{n} = [L{n}; Lz + Lbf + La + Lv, 5];  
                                
                [Lz,~,cp(n)] = update_cp(f,bf,wlnmu,cp(n),msk);
                L{n}         = [L{n}; Lz + Lbf + La + Lv, 2];     
                
                Update = []; ov = [];
            end     
        end    
        
        %------------------------------------------------------------------
        % Update statistics for TPM
        %------------------------------------------------------------------
        
        if do.mu && N>1                                 
            
            % Get likelihoods----------------------------------------------
            
            [~,~,~,Pf] = update_cp(f,bf,wlnmu,cp(n),msk);  
            Pf         = reshape(Pf,[d K]);
            
            % Warp likelihoods to template
            E               = spm_dexpm(a{n},B);
            Affine          = Mmu\E*Mf; 
            invAffine       = inv(Affine);           
            y1              = spm_diffeo('comp',theta,affine_transf(invAffine,identity(d)));
            if d(3)==1, 
                y1(:,:,:,3) = 1; 
            end            
            Pf              = warp(Pf,y1,bs); y1 = []; 
            
            % Scale warped likelihoods by Jacobian determinants
            dt = spm_diffeo('det',Jtheta)*abs(det(invAffine(1:3,1:3)));
            Pf = bsxfun(@times,Pf,dt); dt = [];
            
            % Solve dF/dmu=0 for mu----------------------------------------
            
            % Get logs of mixing weights
            c1 = cp(n).lnw;
            c1 = c1';
            c1 = repmat(c1,1,prod(d)); % [KxNf]
            
            % Get logs of template
            a1 = reshape(lnmu,[prod(d) K])';
            
            ac  = bsxfun(@plus,a1,c1); a1 = [];   
            eac = exp(ac); ac= [];

            % Get likelihoods
            b1   = reshape(Pf,[prod(d) K]); Pf = [];
            b1   = b1';
            eacb = bsxfun(@times,eac,b1); b1= [];

            r1                = bsxfun(@rdivide,eacb,sum(eacb,1)); eacb = [];
            r1(~isfinite(r1)) = 0;   
            
            munum = munum + r1; r1 = []; % Add to sum over subjects

            pn                = 1./sum(eac,1); eac = [];
            cn                = exp(c1); c1 = [];
            p1                = bsxfun(@times,pn,cn); pn = []; cn = [];
            p1(~isfinite(p1)) = 0; 
            
            muden = muden + p1; p1 = []; % Add to sum over subjects
        end
        
        % Save data--------------------------------------------------------
        save_to_nii(n,pthv,v);
        
        if debuglevel>2
            fprintf('cp(%i).w: %s\n',n,sprintf('%3.3f ',exp(cp(n).lnw)));
        end
    end
       
    %----------------------------------------------------------------------
    % Update template (mu)
    %----------------------------------------------------------------------            
    
    if do.mu && N>1     
        
        % Smooth template--------------------------------------------------
        [munum,muden] = smooth_template(munum,muden,d,mufwhm);
        
        % Update template--------------------------------------------------
        lnmu = log(munum./muden + eps); munum = []; muden = []; % eps can be considered a Dirichlet prior
        lnmu = reshape(lnmu',[d K]);
        
        if sum(~isfinite(lnmu(:))), warning('sum(~isfinite(lnmu(:)))'); end                                
        
        % Compute objective------------------------------------------------                
        parfor (n=1:N,runpar)
%         for n=1:N
            [f,v,msk] = load_from_nii(n,pthf,pthv,d,C,ct);
            
            % Re-generate bias field---------------------------------------
            bf = get_bf(chan{n},d);
            
            % Warp temlate to subject--------------------------------------
            if vconv(n), phi = make_deformation(v,prm,int_args,Greens);            
            else         phi = make_deformation(v,prm0,int_args,Greens0); end     
            E                = spm_dexpm(a{n},B);
            Affine           = Mmu\E*Mf;                           
            y1               = affine_transf(Affine,phi);                    
            wlnmu            = warp(lnmu,y1,bs); y1 = [];                                                  
            
            % Objective function-------------------------------------------
            if vconv(n), [Lz,Lbf,La,Lv] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm,msk);            
            else         [Lz,Lbf,La,Lv] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm0,msk); end
            L{n}                        = [L{n}; Lz + Lbf + La + Lv, 6];
                                
            % Update cluster parameters------------------------------------
            [Lz,~,cp(n)] = update_cp(f,bf,wlnmu,cp(n),msk);       
            L{n}         = [L{n}; Lz + Lbf + La + Lv, 2];      
            
            f = []; v = []; bf = []; phi = []; wlnmu = [];
        end
    end                        
    
    %----------------------------------------------------------------------
    % Update intensity prior hyper-parameters
    %----------------------------------------------------------------------    
    
    if do.pr && N>1 
        
        % Update intensity prior hyper-parameters--------------------------
        pr = update_pr(cp);
        for n=1:N
            cp(n).pr = pr;
        end        
        
        % Compute objective------------------------------------------------                
        parfor (n=1:N,runpar)
%         for n=1:N
            [f,v,msk] = load_from_nii(n,pthf,pthv,d,C,ct);
            
            % Re-generate bias field---------------------------------------
            bf = get_bf(chan{n},d);
            
            % Warp temlate to subject--------------------------------------
            if vconv(n), phi = make_deformation(v,prm,int_args,Greens);            
            else         phi = make_deformation(v,prm0,int_args,Greens0); end     
            E                = spm_dexpm(a{n},B);
            Affine           = Mmu\E*Mf;                           
            y1               = affine_transf(Affine,phi);                    
            wlnmu            = warp(lnmu,y1,bs); y1 = [];                                                
            
            % Objective function-------------------------------------------
            if vconv(n), [Lz,Lbf,La,Lv] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm,msk);            
            else         [Lz,Lbf,La,Lv] = obj_fun(do,f,bf,wlnmu,cp(n),chan{n},a{n},R,v,prm0,msk); end
            L{n}                        = [L{n}; Lz + Lbf + La + Lv, 7];
                                
            % Update cluster parameters------------------------------------
            [Lz,~,cp(n)] = update_cp(f,bf,wlnmu,cp(n),msk);       
            L{n}         = [L{n}; Lz + Lbf + La + Lv, 2];    
            
            f = []; v = []; bf = []; phi = []; wlnmu = [];
        end
    end  
    
    %----------------------------------------------------------------------
    % Plot, check stopping criteria, save, etc.
    %----------------------------------------------------------------------
    
    % Display results------------------------------------------------------    
    if debuglevel>1
        visualize_progress(fig,d,K,rndn,lnmu,pthv,chan,zix);
    end
        
    % Save most recent TPM-------------------------------------------------
    if do.writemu
        write_mu(tag,lnmu,Mmu,d,K);
    end
    
    % Plot lower bound-----------------------------------------------------
    if debuglevel
        plot_L(L,itmain,fig);  
    end
    
    % Check convergence----------------------------------------------------
    gL = globalL(gL,L);                  
    if sum(~isfinite(gL(2:end))), warning('~isfinite(gL)'); end

    fprintf('%d %g %g\n',itmain,gL(end),abs((gL(end) - gL(end - 1))/gL(end)));
    
    if abs((gL(end) - gL(end - 1))/gL(end))<tol && itmain>30
        fprintf('Algorithm converged in %d iterations.\n',itmain)                        
        break;
    end
end
%==========================================================================

%==========================================================================
function [pthf,C,N,V] = get_pth2f(imdir,N,C)
folder    = dir(imdir);         % folder with subfolders containing multi-channel data of subjects
folder    = folder(3:end);
dirflag   = [folder.isdir];
subfolder = folder(dirflag);    % subfolders (S1,S2,...)
N1        = numel(subfolder);
if N>N1
    N = N1; % Number of subjects
end

fname = dir(fullfile(imdir,subfolder(1).name,'*.nii'));
C1    = numel(fname);
if C>C1
    C = C1; % Number of channels
end

pthf = cell(N,C);
for n=1:N
    folder = fullfile(imdir,subfolder(n).name);
    fname  = dir(fullfile(folder,'*.nii'));
    for c=1:C
        pthf{n,c} = fullfile(folder,fname(c).name);
    end
end

V = get_all_V(pthf);

fprintf('Loading data from %d subjects having %d channels each\n',N,C); 
%==========================================================================

%==========================================================================
function V = get_all_V(pthf)
[N,C] = size(pthf);

V = cell(N,C);
for c=1:C
    for n=1:N
        V{n,c} = spm_vol(pthf{n,c});
    end
end
%==========================================================================

%==========================================================================
function write_mu(tag,lnmu,Mmu,d,K)
lnmu = reshape(lnmu,[prod(d) K]);
mu   = softmax(lnmu);
mu   = reshape(mu,[d K]);

f     = fullfile('results',tag);
fname = fullfile(f,[tag '__tpm.nii']);
create_nii(fname,mu,Mmu,'float32','TPM');   
%==========================================================================

%=======================================================================        
function [Lz,Lbf,La,Lv] = obj_fun(do,f,bf,lnmu,cp,chan,a,R,v,prm,msk,doLz)
if nargin<12, doLz = true; end

Lz = 0; Lbf = 0; La = 0; Lv = 0;

if doLz
    % E[ln(p(X|Z,µ,Σ,b))] + E[ln(p(Z|π,w,v,a))] + E[ln(p(µ,Σ))] - E[ln(q(Z))] - E[ln(q(µ,Σ))]
    Lz = update_cp(f,bf,lnmu,cp,msk);
end

if do.bf && nargout > 1    
    % ln(p(b))
    Lbf = -0.5*chan(1).T(:)'*chan(1).C*chan(1).T(:);
end

if do.a && nargout > 2
    % ln(p(a))
    La = -0.5*a'*R*a;   
end
if do.v && nargout > 3
    % ln(p(v))
    u  = spm_diffeo('vel2mom',v,prm);
    Lv = -0.5*sum(sum(sum(sum(double(u).*double(v)))));
end     
%=======================================================================

%==========================================================================
function [f,v,msk] = load_from_nii(n,pthf,pthv,d,C,ct)
f   = zeros([prod(d) C],'single');
msk = zeros([prod(d) C],'logical');
for c=1:C
    Nii      = nifti(pthf{n,c});
    f(:,c)   = reshape(Nii.dat(:,:,:),[],1);  
    
    msk(:,c) = get_msk(f(:,c),ct);
end
msk = reshape(msk,[d C]);

v   = nifti(pthv{n});
v   = single(v.dat(:,:,:,:));
%==========================================================================

%==========================================================================
function save_to_nii(n,pthv,v)
Nii              = nifti(pthv{n});
Nii.dat(:,:,:,:) = v;
%==========================================================================

%==========================================================================
function visualize_progress(fig,d,K,rndn,lnmu,pthv,chan,zix)
C = numel(chan{1});

% Velocity fields
set(0,'CurrentFigure',fig{2});
clf(fig{2}); 

cnt = 1;
for f1=1:numel(rndn)
    Nii = nifti(pthv(rndn(f1)));
    v   = Nii.dat(:,:,:,:);

    for f2=1:3
        subplot(numel(rndn),3,cnt)
        imagesc(v(:,:,zix,f2)'); axis image xy off; title(['v' num2str(f2) ', n=' num2str(rndn(f1))]); colormap(gray); colorbar;
        cnt = cnt + 1;
    end
end
drawnow

% Bias fields
set(0,'CurrentFigure',fig{3});
clf(fig{3}); 

cnt = 1;
for f1=1:numel(rndn)
    bf = get_bf(chan{rndn(f1)},d);
    bf = reshape(bf,[d C]);

    c = 1;
    
    bf1 = bf(:,:,floor((size(bf,3) + 1)/2),c);
    bf2 = squeeze(bf(floor((size(bf,1) + 1)/2),:,:,c));
    bf3 = squeeze(bf(:,floor((size(bf,2) + 1)/2),:,c));  

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

% Template
set(0,'CurrentFigure',fig{4});
clf(fig{4}); 

lnmu = reshape(lnmu,[prod(d) K]);
mu   = softmax(lnmu);
mu   = reshape(mu,[d K]);

K1 = floor(sqrt(K));
K2 = ceil(K/K1);      
for k=1:K
    subplot(K1,K2,k);
    imagesc(mu(:,:,zix,k)'); axis image xy off; title(['mu, k=' num2str(k)]); colormap(gray);
end  
drawnow
%==========================================================================

%==========================================================================
function gL = globalL(gL,L,verbose)
if nargin<3, verbose = true; end

sL = zeros(size(L{1},1),1);
for n=1:numel(L)
    sL = sL + L{n}(:,1);
end
    
gL(end + 1) = sL(end);

if gL(end) - gL(end - 1) < 0 && verbose
    fprintf(2,'Global LB decreased by %1.2f\n',gL(end) - gL(end - 1));
end
%==========================================================================

%==========================================================================
function plot_L(L,iter,fig,ix,verbose)
if nargin<4, ix      = 1; end
if nargin<5, verbose = false; end

sL = zeros(size(L{1},1),1);
for n=1:numel(L)
    sL = sL + L{n}(:,1);
end

if any(diff(round(sL,8)) < 0) && verbose
    fprintf(2,'any(diff(sL) < 0)\n');
end

set(0,'CurrentFigure',fig{1});
clf(fig{1}); 

sL     = sL(2:end);
col    = L{1}(2:end,2);
colmap = {'m.','g.','y.','b.','r.','c.','k.'};

hold on
plot(1:numel(sL),sL,'k-','LineWidth',1);
for i=1:numel(sL)
    plot(i,sL(i,ix),colmap{col(i)},'markersize',10);
end
hold off
xlabel('Convergence (global)')    
title(['iter=' num2str(iter)])
drawnow          
%==========================================================================

%==========================================================================
function [fig,rndn,z] = init_fig(N,debuglevel,d,ixstart)
fig = cell(5,1);

ixstart = 1 + 5*(ixstart - 1);

if debuglevel>3
    showwarp = 1;
else
    showwarp = 0;
end

if debuglevel>1
    showres = 1;
else
    showres = 0;
end

if debuglevel
    showL = 1;
else
    showL = 0;
end

if showres && showL && showwarp
    for i=1:5
        fig{i} = figure(ixstart + i - 1); 
        clf(fig{i});     
    end    
elseif showres && ~showL && ~showwarp
    for i=1:5
        if i==1 || i==5
            close(figure(ixstart + i - 1));
        else
            fig{i} = figure(ixstart + i - 1); 
            clf(fig{i});
        end
    end    
elseif showres && showL && ~showwarp
    for i=1:5
        if i==5
            close(figure(ixstart + i - 1));
        else
            fig{i} = figure(ixstart + i - 1); 
            clf(fig{i});           
        end
    end
elseif ~showres && showL && showwarp
    for i=1:5
        if i==1 || i==5
            fig{i} = figure(ixstart + i - 1); 
            clf(fig{i});
        else
           close(figure(ixstart + i - 1));
        end
    end
elseif ~showres && ~showL && showwarp
    for i=1:5
        if i==5
            fig{i} = figure(ixstart + i - 1); 
            clf(fig{i});
        else
           close(figure(ixstart + i - 1));
        end
    end    
elseif ~showres && showL && ~showwarp
    for i=1:5
        if i==1
            fig{i} = figure(ixstart + i - 1); 
            clf(fig{i});
        else
           close(figure(ixstart + i - 1));
        end
    end     
else
    for i=1:5
        close(figure(ixstart + i - 1));
    end
end
drawnow;

rndn = datasample(1:N,min(N,4),'Replace',false);

z = floor((d(3) + 1)/2);
%==========================================================================

%==========================================================================
function [L,gL] = init_L(N)
L = cell(N,1);
for n=1:N
    L{n} = [-Inf 1];
end
gL = -Inf;
%==========================================================================

%==========================================================================
function show_mu(fig,K,lnmu,r,d,zix)
set(0,'CurrentFigure',fig);       
                      
lnmu = reshape(lnmu,[prod(d) K]);
mu   = softmax(lnmu);
mu   = reshape(mu,[d K]);

for k=1:K
    subplot(2,K,k  ); imagesc(mu(:,:,zix,k)'); axis image xy off; colormap(gray); title(['mu' num2str(k)]);
    subplot(2,K,k+K); imagesc(r(:,:,zix,k)');  axis image xy off; colormap(gray); title(['r' num2str(k)]);
end
drawnow
%==========================================================================