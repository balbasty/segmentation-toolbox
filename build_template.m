function build_template

% PRIO
% Add back missing 
%  --problem with convergence of bf (when lkp>1?), probably due to
%    decreasing ll...
%  --include missing ll for ml
%
% -Clean up repo + add readme
% -Labelled resps
% -Improved registration: diffeo + reparameterise + average correct
%
% TODO
% -Add denoising
% -Add superres
% -Add remove corruped
% -Add DICOM convert
% -Introduce aux
% -Use K - 1 classes for template?
% -Investigate shoot_template prm
% -Create pars.seg.
% -K.rem,K.keep,K.lkp
%
% CLUSTER
% -Reintroduce split?
% -Same results, cluster vs non-cluster
%

%--------------------------------------------------------------------------
addpath(genpath('./code'))
addpath('/cherhome/mbrud/dev/distributed-computing/')
addpath('/cherhome/mbrud/dev/auxiliary-functions/')

%--------------------------------------------------------------------------
NAME    = 'CT-lesion-local';
TESTING = 0;

%-------------------------------------------------------------------------- 
im = {};
         
if TESTING
    if     TESTING==1 || TESTING==4, S = 1;
    elseif TESTING>1,                S = 8; 
    end
    
    K           = 6;               
    im{end + 1} = {'/data/mbrud/images/IXI-subjects/',...
                   S,'MRI',[],4,'mean',''};   
%     im{end + 1} = {'/data/mbrud/images/OASIS-long-noneck/',...
%                    S,'MRI',[],4,'mean',''};                  
%     im{end + 1} = {'/data/mbrud/images/CT-healthy-noneck-den/',...
%                    S,'CT',[],4,'mean',''};                 
else
    K           = 8;
    im{end + 1} = {'/data/mbrud/images/CT-CHROMIS-noneck-den/',...
                   Inf,'CT',[],2,'mean',''};    
    im{end + 1} = {'/data/mbrud/images/CT-healthy-noneck-den/',...
                   Inf,'CT',[],2,'mean',''};                        
%     im{end + 1} = {'/data/mbrud/images/CT-big-lesions/',...
%                    Inf,'CT',[],3,'total',''};
%     im{end + 1} = {'/data/mbrud/images/OASIS-long-noneck/',...
%                    30,'MRI',[],3,'total',''};                  
end

%--------------------------------------------------------------------------
if TESTING==1 || TESTING==2
    NAME = 'test-local';
elseif TESTING==3
    NAME = 'test-holly';
elseif TESTING==4
    NAME = 'test-preproc';    
end

dir_output                = '/data-scratch/mbrud/data/';
dir_template              = '/data/mbrud/templates/';
[dir_output,dir_template] = append_dir(dir_output,dir_template,NAME);

%--------------------------------------------------------------------------
holly = struct;

holly.server.ip     = 'holly';
holly.server.login  = 'mbrud';
holly.server.folder = fullfile('/scratch',dir_output(15:end),'cluster');
holly.client.folder = fullfile(dir_output,'cluster');

if TESTING==1 || TESTING==2 || TESTING==4
    holly.server.ip  = '';   
    
    if TESTING==1 || TESTING==4, 
        holly.client.workers = 0;
    elseif TESTING==2
        holly.client.workers = Inf;
    end
end

holly.server.ip      = '';   
holly.client.workers = Inf;

holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = '/home/mbrud/dev/build-template';

holly.translate     = {'/data-scratch/mbrud/' '/scratch/mbrud/'};
holly.restrict      = 'char';
holly.clean         = true;
holly.verbose       = false;

holly.job.est_mem   = true;
holly.job.batch     = true;
holly.job.mem       = '6G';
holly.job.use_dummy = true;

%--------------------------------------------------------------------------
pars.do_segment      = true;
pars.do_preproc      = false;
pars.do_ml           = false;
pars.do_bf           = true;
pars.do_def          = true;
pars.do_wp           = true;
pars.do_write_res    = false;
pars.do_push_resp    = false;
pars.do_missing_data = false;
pars.do_old_segment  = false;

if TESTING==4
    pars.do_preproc = true;
    pars.do_segment = false;
end

pars.preproc.reg_and_reslice = true;
pars.preproc.realign2mni     = true;
pars.preproc.crop            = true;
pars.preproc.rem_neck        = true;

pars.kmeans_dist     = 'cityblock';
pars.niter           = 50;
pars.dir_output      = dir_output;

pars.print_ll        = false;
pars.print_seg       = false;    

do_shrink           = true;
tol                 = 1e-4;
do_avg_template_dim = true;

do_show_seg     = false;   
do_show_results = 3;
if TESTING==1
    do_show_seg    = true;
    pars.print_ll  = true;
    pars.print_seg = true;
elseif TESTING==2
    do_show_seg    = false;
    pars.print_ll  = false;
    pars.print_seg = false;    
end

%--------------------------------------------------------------------------
M               = numel(im);
V               = cell(1,M);
for m=1:M, V{m} = get_V(im{m}); end

%--------------------------------------------------------------------------
pth_template = '';
% pth_template = '/data-scratch/mbrud/data/build-template-ct/template/template.nii';

% Uncomment below to use predefined templates
% pth_template = '/home/mbrud/Dropbox/PhD/Data/log-template/logTPM.nii';
% pars.lkp     = [1 1 2 2 3 3 4 4 5 5 5 6 6];
% pth_template = fullfile(get_pth_dropbox,'/PhD/Data/CB-TPM/BlaiottaTPM.nii');
% pars.lkp     = [1 1 2 2 3 3 4 4 5 5 6 6 7 7];

[pth_template,uniform,dt] = init_template(pth_template,V,K,dir_template,do_avg_template_dim); 

%--------------------------------------------------------------------------
pth_prior = '';
% pth_prior = '/data-scratch/mbrud/data/build-template-ct/template/prior-CT-big-lesions.mat';

%------------------------------------------------------
[obj,niter,fig,rand_subjs] = init_obj(V,im,pth_template,pth_prior,uniform,do_show_seg,do_show_results,pars,dt,TESTING);

%-------------------------------------------------------------------------- 
holly = distribute_default(holly);

%--------------------------------------------------------------------------
print_algorithm_progress('started');

L = -Inf;
for iter=1:niter        
    
    if niter>1
        % Some parameters of the obj struct are changed depending on iteration 
        % number (only for building templates)
        %----------------------------------------------------------------------    
        obj = modify_obj(obj,iter);
    end       
    
    % Segment a bunch of subjects 
    %----------------------------------------------------------------------
    [obj,ix]    = unfold_cell(obj,2);
    [holly,obj] = distribute(holly,'update_subject','inplace',obj,pth_template,fig);
    obj         = fold_cell(obj,ix);

    % Check if any subjects have status~=0
    %----------------------------------------------------------------------
    print_jobs_failed(obj);
    
    if niter>1
        % Update template
        %------------------------------------------------------------------
        L = update_template(L,pth_template,obj,iter);                              
    end        
    
    if niter>1 && do_shrink
        % Automatically decrease the size of the template based on bounding
        % boxes computed for each pushed responsibility image
        %------------------------------------------------------------------
        shrink_template(pth_template,obj,iter);    
    end
    
    if niter>1
        % Update Gaussian-Wishart hyper-parameters
        %------------------------------------------------------------------
        obj = update_intensity_prior(obj,dir_template,iter);
    end
       
    % Save object file
    %----------------------------------------------------------------------
    save(fullfile(dir_template,'obj.mat'),'obj')
    
    if niter>1
        % Some verbose
        %----------------------------------------------------------------------
        if do_show_results>0, plot_ll(fig{5},L); end
        if do_show_results>1, show_template(fig{6},pth_template); end
        if do_show_results>2, show_resp(fig{7},obj,rand_subjs); end      

        tol1 = abs((L(end - 1)*(1 + 10*eps) - L(end))/L(end));    
        fprintf('%2d | L = %0.0f | diff = %0.7f | tol = %0.5f\n',iter,L(end),L(end) - L(end - 1),tol1);  
        
        if tol1<tol
           break 
        end
    end
end

print_algorithm_progress('finished',iter);
%==========================================================================

%==========================================================================
function obj = modify_obj(obj,iter)
M = numel(obj);    
for m=1:M
    S         = numel(obj{m});    
    sum_bf_dc = 0;
    cnt_S     = 0;
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==1 
            obj{m}{s}.do_def       = false;
            obj{m}{s}.do_bf        = false;
            obj{m}{s}.do_push_resp = true;
            obj{m}{s}.do_write_res = false;
            
            obj{m}{s}.niter  = 1;                 
            obj{m}{s}.nsubit = 1;
            obj{m}{s}.nitgmm = 1;
        end

        if iter==2
            obj{m}{s}.nsubit  = 8;
            obj{m}{s}.nitgmm  = 20;    
            if strcmp(obj{m}{s}.modality,'MRI')
                obj{m}{s}.do_bf = true;            
            end
            obj{m}{s}.uniform = false;
        end

        if iter>=3
            obj{m}{s}.do_def = true;         
            obj{m}{s}.reg    = obj{m}{s}.reg0;            
            scal             = 2^max(11 - iter,0);       
            %prm([5 7 8]) = param([5 7 8])*scal;
            obj{m}{s}.reg(3) = obj{m}{s}.reg(3)*scal;
        end
        
        % Sum bias field DC components
        sum_bf_dc = sum_bf_dc + obj{m}{s}.bf_dc;
        cnt_S     = cnt_S + 1;
    end
    
    % Average of bias field DC components
    avg_bf_dc = sum_bf_dc/cnt_S;
    
    % Set average bias field DC component
    for s=1:S 
        obj{m}{s}.avg_bf_dc = avg_bf_dc; 
    end
end
%==========================================================================

%==========================================================================
function shrink_template(pth_template,obj,iter,verbose)
if nargin<4, verbose = true; end

M   = numel(obj);
bb  = [];
cnt = 1;
for m=1:M
    S = numel(obj{m});    
    for s=1:S            
        if obj{m}{s}.status==0                
            bb(:,:,cnt) = obj{m}{s}.bb_push;
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
    subvol(V0(k),nbb','tmp');        
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
function show_template(fig,pth_template)
set(0,'CurrentFigure',fig);     

Nii = nifti(pth_template);
b   = exp(Nii.dat(:,:,:,:));
b   = bsxfun(@rdivide,b,sum(b,4));
d   = size(b);
K   = d(4);                                  
for k=1:K         
    subplot(3,K,k);
    slice = b(:,:,floor(d(3)/2) + 1,k);
    imagesc(slice'); axis image xy off; title(['k=' num2str(k)]); colormap(gray);               
   
    subplot(3,K,K + k);
    slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
    imagesc(slice); axis image xy off; colormap(gray);   

    subplot(3,K,2*K + k);
    slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
    imagesc(slice'); axis image xy off; colormap(gray);   
end 
drawnow
%==========================================================================

%==========================================================================
function show_resp(fig,obj,rand_subjs)
M = numel(obj);
set(0,'CurrentFigure',fig);       

cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        K = numel(obj{m}{s}.pth_resp);
        
        for k=1:K    
            Nii = nifti(obj{m}{s}.pth_resp{k});
            img = Nii.dat(:,:,:);
            dm  = size(img);
            zix = floor(dm(3)/2) + 1;
            img = img(:,:,zix);
        
            subplot(M*numel(rand_subjs{1}),K,cnt);
            imagesc(img'); axis image xy off; colormap(gray);
            title(['q_{' num2str(m), ',' num2str(s) ',' num2str(k) '}']);
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
function [dir_output,dir_template] = append_dir(dir_output,dir_template,name)
pth        = fileparts(dir_output);
dir_output = fullfile(pth,['build-template-' name]);

pth          = fileparts(dir_template);
dir_template = fullfile(pth,name);
%==========================================================================