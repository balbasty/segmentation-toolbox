function test_GaussiansFromSuffStats

addpath(genpath('../code'));

% Get paths
S           = Inf;
im          = {};
im{end + 1} = {'/home/mbrud/Dropbox/PhD/Data/2D-Data/TEST-MOG/IXI',S,'MRI','healthy',''};
im{end + 1} = {'/home/mbrud/Dropbox/PhD/Data/2D-Data/TEST-MOG/OASIS',S,'MRI','healthy',''};
im{end + 1} = {'/home/mbrud/Dropbox/PhD/Data/2D-Data/TEST-MOG/CHROMIS',S,'CT','healthy',''};

% Parameters
obj.preproc.do_preproc    = false;
obj.preproc.rem_corrupted = false;
obj.num_workers           = false;
obj.run_on_holly          = false;

obj.tol   = 1e-4;
obj.niter = 500;

obj.missing_data = false;
obj.dowp         = false; % false for uniform TPM?
obj.wp_reg       = 1;

iterstop = 10;
mrange   = [1];
Krange   = [6];
vbrange  = [true];

verbose = true;

ok = [];
for m=mrange
    
    V = load_and_process_images(obj,im(m));
    
    for K=Krange         

        for vb=vbrange
            fprintf('\n= m=%d, K=%d, vb=%d ===============================\n',m,K,vb);

            % Create buffer
            lkp          = 1:K;
            missing_data = obj.missing_data;
            wp_reg       = obj.wp_reg;
            dowp         = obj.dowp;
            tol          = obj.tol;
            niter        = obj.niter;

            cast = @uint8;
            typ  = 'uint8';
            buf  = struct;

            N  = numel(V{1}{1});
            nz = 1;

            s0 = zeros(1,N);
            s1 = zeros(1,N);
            S2 = zeros(1,N);
            Nm   = 0;
            for z=1:nz       

                for n=1:N
                    Nii = nifti(V{1}{1}(n).fname);
                    f   = squeeze(single(Nii.dat(:,:,:)));

                    buf(z).msk{n} = get_msk(f,im{m}{3});
                    buf(z).nm{n}  = nnz(buf(z).msk{n});
                end

                if ~missing_data
                    msk = true;
                    for n=1:N
                        msk = msk & buf(z).msk{n};
                    end
                    
                    for n=1:N
                        buf(z).msk{n} = msk;
                        buf(z).nm{n}  = nnz(buf(z).msk{n});
                    end
                end
                clear msk
                
                for n=1:N
                    Nii = nifti(V{1}{1}(n).fname);
                    f   = squeeze(single(Nii.dat(:,:,:)));

                    buf(z).f{n}  = f(buf(z).msk{n});  
                    buf(z).bf{n} = ones(size(buf(z).f{n}),'single');  

                    s0(n) = s0(n) + buf(z).nm{n};
                    s1(n) = s1(n) + sum(buf(z).f{n});
                    S2(n) = S2(n) + sum(buf(z).f{n}.^2);
                end    

                code = zeros([numel(buf(z).msk{1}) 1],typ);
                for n=1:N
                    code = bitor(code,bitshift(feval(cast,buf(z).msk{n}(:)),(n-1)));
                end
                buf(z).code = code;

                buf(z).Nm = nnz(code~=0);
                Nm        = Nm + buf(z).Nm;        

                buf(z).dat = ones([buf(z).Nm K],'single');
            end

            %% Run algorithm
               
            mog     = [];
            mog.vb  = vb;
            mog.vr0 = diag(S2./s0 - (s1./s0).^2)/K^2;
            clear s0 s1 S2
            
            wp = ones(1,K)/K;
            mg = ones(1,K);

            L  = -Inf;
            lb = -Inf;
            lh = -Inf;
            lq = -Inf(1,3);
            Lv = [];

            t = tic;
            for i=1:niter
                oL  = L;   
                olb = lb;
                olh = lh;
                olq = lq;
                L   = 0;

                % M-step   
                %----------------------------------------------------------------------

                if i==1                               
                    mom       = spm_kmeans2mom(buf,K,verbose);
                    [~,~,mog] = spm_GaussiansFromSuffStats(mom,mog);  
                else        
                    [~,~,mog,s0] = spm_GaussiansFromSuffStats(mom,mog);                  

                    if dowp
                        % Update mixing weights
                        mgm = 0;
                        for z=1:nz       
                            tpm = double(buf(z).dat);
                            s   = 1./(tpm*wp');                                        
                            mgm = mgm + s'*tpm; 
                        end
                        for k1=1:K
                            wp(k1) = (sum(s0(lkp==k1)) + wp_reg*1)/(mgm(k1) + wp_reg*K);
                        end
                        wp = wp/sum(wp);    
                    end
                end    

                % E-step        
                %----------------------------------------------------------------------

                lq  = 0;
                mom = mom_struct(K,N); 
                for z=1:nz        
                    cr = zeros(numel(buf.msk{1}),N);
                    for n=1:N
                        cr(buf.msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n});
                    end

                    B = double(buf(z).dat);
                    
                    % Calculate responsibilities
                    [Q,dlq] = latent(buf(z),B,mg,mog,wp,lkp,cr);
                    lq      = lq + dlq;

                    % Calculate sufficient statistics                
                    mom = spm_SuffStats(cr,Q,buf(z).code,mom);
                end    
                L = L + sum(lq);    

                [lb,lh] = spm_GaussiansFromSuffStats(mom,mog);  
                L       = L + lb + lh;
                Lv      = [Lv,L];

                my_fprintf(verbose,i,wp,mog,lb,lh,lq,olb,olh,olq,L,oL);

                if i>iterstop && L-oL<tol*Nm       
                    break;
                end  
            end
            t = toc(t);

            my_fprintf(true,i,wp,mog,lb,lh,lq,olb,olh,olq,L,oL,t);

            ok = [ok,min(diff(single(Lv)))];

            if verbose
                % Display estimated responsibilities    
                figure(668);
                d   = [size(buf(1).msk{1}) 1];
                zix = floor(d(3)/2) + 1;
                Q   = reshape(Q,[d K]);                   
                K1  = floor(sqrt(K));
                K2  = ceil(K/K1); 
                for k=1:K
                  subplot(K1,K2,k);
                  imagesc(Q(:,:,zix,k)'); axis image xy off; colormap(gray);
                end
            end
        end
    end
end

fprintf('\n==================================\n')
disp(num2str(ok)); % should be greater than zero
fprintf('==================================\n')
%=======================================================================

%======================================================================= 
function my_fprintf(verbose,i,wp,mog,lb,lh,lq,olb,olh,olq,L,oL,t)
if nargin<13, t = 0; end
if verbose
    % Cast to single data type precision
    lb  = single(lb);
    lh  = single(lh);
    lq  = single(lq);
    olb = single(olb);
    olh = single(olh);
    olq = single(olq);
    L   = single(L);
    oL  = single(oL);
    
    fprintf('= iter=%d ================================\n',i);
    % Print mixing proportions
    fprintf('wp    = [%0.6f, %s%0.6f]\n',wp(1),sprintf('%0.6f, ',wp(2:end-1)),wp(end)); 
    if isfield(mog,'pr')
        % Print posterior
        b  = mog.po.b;
        n  = mog.po.n;
        m  = mog.po.m;
        N  = size(m,1);    

        fprintf('po.b  = [%0.2f, %s%0.2f]\n',b(1),sprintf('%0.2f, ',b(2:end-1)),b(end)); 
        fprintf('po.n  = [%0.2f, %s%0.2f]\n',n(1),sprintf('%0.2f, ',n(2:end-1)),n(end));     
        for n=1:N
           fprintf('po.m%d = [%0.2f, %s%0.2f]\n',n,m(n,1),sprintf('%0.2f, ',m(n,2:end-1)),m(n,end)); 
        end        
        fprintf('lb    = %0.2f,\t(lb - olb = %0.2f)\n',lb,lb - olb);    
        fprintf('lh    = %0.2f,\t(lh - olh = %0.2f)\n',lh,lh - olh);                
        fprintf('lq    = [%0.2f, %0.2f, %0.2f],\t(lq - olq = [%0.2f, %0.2f, %0.2f])\n',lq(1),lq(2),lq(3),lq(1) - olq(1),lq(2) - olq(2),lq(3) - olq(3));     
    else
        % Print mean
        mn = mog.mn;
        N  = size(mn,1);
        for n=1:N
           fprintf('mn%d   = [%0.2f, %s%0.2f]\n',n,mn(n,1),sprintf('%0.2f, ',mn(n,2:end-1)),mn(n,end)); 
        end
    end
    % Print lower bound
    diff1 = L - oL;
    if diff1>=0
        fprintf('L     = %0.2f,\t(L - oL = %0.7f)\n',L,diff1);
    else
        fprintf(2,'L     = %0.2f,\t(L - oL = %0.7f)\n',L,diff1);
    end
    if t
        % Print elapsed time
        fprintf('t     = %0.2f s\n',t);
    end
end
%======================================================================= 