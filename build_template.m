function build_template

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
dir_output = '/home/smajjk/Data/tmp-build-tpm';

pth2_distributed_toolbox = '../distributed-computing';
pth2_auxiliary_functions = '../auxiliary-functions';

holly_server_login   = 'mbrud';
holly_matlab_add_src = '/home/mbrud/dev/segmentation-toolbox';
holly_matlab_add_aux = '/home/mbrud/dev/auxiliary-functions';

% addpath
%--------------------------------------------------------------------------

addpath(genpath('./code'))
addpath(pth2_distributed_toolbox)
addpath(pth2_auxiliary_functions)

test_level = 2; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 8 subjects (holly)
m = 0;

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------

pars.name       = 'CHROMIS';
pars.dir_output = dir_output;
pars.dat        = {};

% Basic test
% -----------------
m = m + 1;
% pars.dat{m}.dir_data = '/home/smajjk/Dropbox/PhD/Data/IXI-test/2d_IXI-T1T2PD_preproc-ra-cr-rn-reg-res-vx';
pars.dat{m}.dir_data = '/home/smajjk/Dropbox/PhD/Data/IXI-test/IXI-T1T2PD_preproc-ra-cr-rn-reg-res-vx';
pars.dat{m}.S = 4;
pars.dat{m}.segment.samp = 1;

% Define a log template
%-----------------
% pars.pth_template = '/mnt/cifs_share/share_data/log-TPMs/CB/logBlaiottaTPM.nii'; % In Ashburner_group shared
% pars.pth_template = '/mnt/cifs_share/share_data/log-TPMs/SPM/logTPM.nii'; % In Ashburner_group shared

pars = pars_default(pars,test_level); 

%--------------------------------------------------------------------------
% Set distribute package parameters
%--------------------------------------------------------------------------

holly               = struct;
holly.server.ip     = 'holly';
holly.server.login  = holly_server_login;
holly.client.folder = fullfile(pars.dir_output,'cluster');
holly.server.folder = holly.client.folder;
holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = holly_matlab_add_src;
holly.matlab.add    = holly_matlab_add_aux;
holly.restrict      = 'char';
holly.clean         = false;
holly.clean_init    = true;
holly.verbose       = false;
holly.job.mem       = '10G';
holly.job.use_dummy = true;

if     test_level==1, holly.server.ip  = ''; holly.client.workers = 0;
elseif test_level==2, holly.server.ip  = ''; holly.client.workers = Inf;
end
% holly.server.ip = ''; holly.client.workers = Inf;
% holly.server.ip = ''; holly.client.workers = 0;

holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Initialise algorithm
%--------------------------------------------------------------------------

pars       = read_images(pars); 
pars       = init_template(pars); 
[obj,pars] = init_obj(pars);

obj = kmeans_on_hist(obj,pars);

%--------------------------------------------------------------------------
% Start the algorithm
%--------------------------------------------------------------------------

print_algorithm_progress('started');

L = -Inf;
for iter=1:pars.niter        
    
    if pars.niter>1
        % Some parameters of the obj struct are changed depending on iteration 
        % number (only for building templates)
        %----------------------------------------------------------------------    
        obj = modify_obj(obj,iter,pars.niter);
    end       
    
    % Segment a bunch of subjects 
    %----------------------------------------------------------------------
    [obj,ix]    = unfold_cell(obj,2);
    [holly,obj] = distribute(holly,'process_subject','inplace',obj,pars.fig);
    obj         = fold_cell(obj,ix);
    
    % Check if any subjects have status~=0
    %----------------------------------------------------------------------
    print_jobs_failed(obj);
    
    if pars.niter>1
        % Update template
        %------------------------------------------------------------------
        L = update_template(L,obj,pars,iter);                              
    end        
    
    if pars.niter>1 && pars.mrf
        % Use a MRF cleanup procedure
        %------------------------------------------------------------------
        mrf_clean_up(pars.pth_template);
    end
        
    if pars.niter>1
        % Automatically decrease the size of the template based on bounding
        % boxes computed for each pushed responsibility image
        %------------------------------------------------------------------
        shrink_template(obj,iter);    
    end
    
    if pars.niter>1 && pars.crop_template==iter
        % Crop template to size of the default SPM template
        %------------------------------------------------------------------
        crop_template(pars.pth_template,iter);
    end
    
    if pars.niter>1
        % Update Gaussian-Wishart hyper-parameters
        %------------------------------------------------------------------
        obj = update_intensity_prior(obj,iter);
    end
       
    % Save obj structs
    %------------------------------------------------------------------
    save(fullfile(pars.dir_template,'obj.mat'),'obj');
    
    if pars.niter>1
        % Some verbose
        %------------------------------------------------------------------
        if pars.verbose>0, plot_ll(pars.fig{5},L); end
        if pars.verbose>1, show_template(pars.pth_template,pars.fig{6}); end
        if pars.verbose>2, show_resp(pars.fig{7},obj,pars); end      
        if pars.verbose>3, show_def(pars.fig{8},obj,pars); end 

        d = abs((L(end - 1)*(1 + 10*eps) - L(end))/L(end));    
        fprintf('%2d | L = %0.0f | d = %0.5f\n',iter,L(end),d);  
        
%         if d<pars.tol
%            break 
%         end
    end
end

print_algorithm_progress('finished',iter);
%==========================================================================

%==========================================================================
function obj = modify_obj(obj,iter,niter)
M = numel(obj);    
for m=1:M
    S         = numel(obj{m});    
    sum_bf_dc = 0;
    cnt_S     = 0;
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==1             
            obj{m}{s}.push_resp.do_push_resp = true;
            obj{m}{s}.write_res.do_write_res = false;
            
            obj{m}{s}.segment.do_def = false;
            obj{m}{s}.segment.do_bf  = false;
            obj{m}{s}.segment.do_wp  = false;
            obj{m}{s}.segment.niter  = 1;                 
            obj{m}{s}.segment.nsubit = 1;
            obj{m}{s}.segment.nitgmm = 1;
        end

        if iter==2            
            obj{m}{s}.uniform = false;  
            
            obj{m}{s}.segment.nsubit = 8;
            obj{m}{s}.segment.nitgmm = 20;  
            obj{m}{s}.segment.do_bf  = obj{m}{s}.segment.do_bf0;                                                  
            obj{m}{s}.segment.do_wp  = obj{m}{s}.segment.do_wp0;    
        end

        if iter>=2
            reg0  = obj{m}{s}.segment.reg0;   
            sched = 2.^fliplr(repelem(0:10,2));
            scal  = sched(min(iter,numel(sched)));   
            
            obj{m}{s}.segment.reg([2 3 4 5]) = reg0([2 3 4 5])*scal;                     
        end
        
        % Sum bias field DC components
        sum_bf_dc = sum_bf_dc + obj{m}{s}.segment.bf_dc;
        cnt_S     = cnt_S + 1;
        
        if iter==niter && obj{m}{s}.image(1).dim(3)>1
            obj{m}{s}.write_res.do_write_res = true;
            obj{m}{s}.write_res.mrf = 2;
            obj{m}{s}.write_res.write_tc(:,[1 2 4]) = true;            
        end
    end
    
    % Average of bias field DC components
    avg_bf_dc = sum_bf_dc/cnt_S;
    
    % Set average bias field DC component
    for s=1:S 
        obj{m}{s}.segment.avg_bf_dc = avg_bf_dc; 
    end
end
%==========================================================================

%==========================================================================
function mrf_clean_up(pth_template,verbose)
if nargin<2, verbose = false; end

Nii = nifti(pth_template);
mat = Nii.mat;
vx  = 1./single(sum(mat(1:3,1:3).^2));
Q   = single(Nii.dat(:,:,:,:));
dm  = size(Q);
Kb  = dm(4);
zix = floor(dm(3)/2) + 1;

% softmax
Q = exp(Q);
Q = bsxfun(@rdivide,Q,sum(Q,4));

nmrf_its = 10;
T        = 2;
G        = T*ones([Kb,1],'single');
P        = zeros(dm,'uint8');

if verbose
    figure(666);
    for k=1:Kb
       subplot(2,Kb,k) 
       imagesc(Q(:,:,zix,k)); axis off image xy; colormap(gray);
    end
end

for iter=1:nmrf_its
    spm_mrf(P,Q,G,vx);
end

P = double(P)/255;

if verbose
    for k=1:Kb
       subplot(2,Kb,Kb + k) 
       imagesc(P(:,:,zix,k)); axis off image xy; colormap(gray);
    end
    drawnow
end

Nii.dat(:,:,:,:) = log(max(P,eps('single')));
%==========================================================================

%==========================================================================
function crop_template(pth_template,iter,verbose)
if nargin<3, verbose = true; end

pth0 = fullfile(spm('dir'),'tpm','TPM.nii');  
V0   = spm_vol(pth0);   
d0   = V0(1).dim;
vx0  = vxsize(V0(1).mat);

V1  = spm_vol(pth_template);
K1  = numel(V1);
d1  = V1(1).dim;
vx1 = vxsize(V1(1).mat);

msk     = d1<d0;
d0(msk) = d1(msk);

sk0 = d0;
sk1 = d1;

bb1 = floor((sk1 - sk0)/2);
bb2 = bb1 + sk0;
bb  = [bb1' bb2'];

for k=1:K1
    spm_impreproc('subvol',V1(k),bb','');        
end

if verbose
    fprintf('%2d | size(otpm) = [%d %d %d] | size(ntpm) = [%d %d %d]\n',iter,d1(1),d1(2),d1(3),d0(1),d0(2),d0(3));
end
%==========================================================================

%==========================================================================
function shrink_template(obj,iter,verbose)
if nargin<3, verbose = true; end

pth_template = obj{1}{1}.pth_template;

M   = numel(obj);
bb  = [];
cnt = 1;
for m=1:M
    S = numel(obj{m});    
    for s=1:S            
        if obj{m}{s}.status==0                
            bb(:,:,cnt) = obj{m}{s}.push_resp.bb;
            cnt         = cnt + 1;
        end
    end
end

mn_bb = min(bb,[],3);
mx_bb = max(bb,[],3);
nbb   = [mn_bb(:,1) mx_bb(:,2)];

V0 = spm_vol(pth_template);
od = V0(1).dim;

for k=1:numel(V0)
    spm_impreproc('subvol',V0(k),nbb','tmp');        
end

delete(pth_template);
[pth,nam,ext] = fileparts(V0(1).fname);
fname         = fullfile(pth,['tmp' nam ext]);
movefile(fname,pth_template);

V  = spm_vol(pth_template);
nd = V(1).dim;   

if verbose
    fprintf('%2d | size(otpm) = [%d %d %d] | size(ntpm) = [%d %d %d]\n',iter,od(1),od(2),od(3),nd(1),nd(2),nd(3));
end
%==========================================================================

%==========================================================================
function print_jobs_failed(obj)
cnt = 0;
for m=1:numel(obj)
    S = numel(obj{m}); 
    for s=1:S
        cnt      = cnt + obj{m}{s}.status;
    end    
end

if cnt
   fprintf('% i job(s) failed!\n',cnt);
end
%==========================================================================

%==========================================================================
function show_resp(fig,obj,pars)
set(0,'CurrentFigure',fig);       

rand_subjs = pars.rand_subjs;

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        K = numel(obj{m}{s}.pth_resp);
        
        for k=1:K    
            Nii = nifti(obj{m}{s}.pth_resp{k});
            img = Nii.dat(:,:,:);
            dm  = size(img);
            if numel(dm)==3
                zix = floor(dm(3)/2) + 1;
                img = img(:,:,zix);            
            end
            
            wp = round(obj{m}{s}.segment.wp(k),3);
            
            subplot(M*numel(rand_subjs{1}),K,cnt);
            imagesc(img'); axis off image xy; colormap(gray);
            title(['q_{' num2str(m), ',' num2str(s) ',' num2str(k) '},w=' num2str(wp)]);
            cnt = cnt + 1;
        end 
    end
end                                  
drawnow
%==========================================================================

%==========================================================================
function show_def(fig,obj,pars)
set(0,'CurrentFigure',fig);       

rand_subjs = pars.rand_subjs;

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        Nii = nifti(obj{m}{s}.pth_vel);
        img = Nii.dat(:,:,:,:);
        dm  = size(img);
        zix = floor(dm(3)/2) + 1; 
        
        for i=1:3
            subplot(M*numel(rand_subjs{1}),3,cnt)
            imagesc(img(:,:,zix,i)'); axis off image xy; colormap(gray); colorbar
            title(['def_{' num2str(m), ',' num2str(s), ',' num2str(i) '}']);
            cnt = cnt + 1;
        end 
    end
end                                  
drawnow
%==========================================================================

%==========================================================================
function print_algorithm_progress(status,iter)
if nargin<2, iter = 0; end

date    = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=================================================\n')  
if iter
    fprintf('| %s | Algorithm %s in %i iterations\n',date,status,iter)
else
    fprintf('| %s | Algorithm %s\n',date,status)
end
fprintf('=================================================\n\n')   
%==========================================================================

%==========================================================================
function plot_ll(fig,L)
set(0,'CurrentFigure',fig);                
plot(0:numel(L(3:end)) - 1,L(3:end),'b-','LineWidth',1);   hold on            
plot(0:numel(L(3:end)) - 1,L(3:end),'b.','markersize',10); hold off  
title('ll')
%==========================================================================

%==========================================================================
function obj = kmeans_on_hist(obj,pars)
K = pars.K;

x = -2000:2000;
h = zeros(1,numel(x));
for m=1:numel(obj)
    if pars.dat{m}.segment.kmeans_hist
        S = numel(obj{m});
        for s=1:S
            img = single(obj{m}{s}.image(1).private.dat(:,:,:));
            msk = msk_modality(img,obj{m}{s}.modality,obj{m}{s}.trunc_ct);
            img = img(msk(:));
            
            h1 = hist(img(:),x);
            h  = h + h1;
        end
    end
end
h = h(900:end)';
x = x(900:end)';

[mg,mn,vr] = spm_imbasics('fit_gmm2hist',h,x,K);

for m=1:numel(obj)
    if pars.dat{m}.segment.kmeans_hist
        S = numel(obj{m});
        for s=1:S
            for k=1:K
                if obj{m}{s}.segment.do_ml                
                    obj{m}{s}.segment.gmm.mn(:,k) = mn(k);
                    obj{m}{s}.segment.gmm.vr(:,:,k) = vr(k);
                else
                    obj{m}{s}.segment.gmm.pr.m(:,k) = mn(k);
                    obj{m}{s}.segment.gmm.pr.b(k) = mean(mg);
                    obj{m}{s}.segment.gmm.pr.n(k) = mean(mg);
                    obj{m}{s}.segment.gmm.pr.W(:,:,k) = 1./(vr(k)*mean(mg));

                    obj{m}{s}.segment.gmm.po.m(:,k) = mn(k);
                    obj{m}{s}.segment.gmm.po.b(k) = mean(mg);
                    obj{m}{s}.segment.gmm.po.n(k) = mean(mg);
                    obj{m}{s}.segment.gmm.po.W(:,:,k) = 1./(vr(k)*mean(mg));
                end      
            end
        end
    end
end
%==========================================================================