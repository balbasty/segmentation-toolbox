function build_template

% 1. Parameterise in log-space
%     --Set tiny and deg to zero for template construction?
% 2. Add back missing 
%     --problem with convergence of reg and bf
%     --include missing ll for ml
%     --artefact: BECAUSE THERE ARE NEG VALUES INTRODUCED WHEN CREATING PREPROC DATA
%
% -Clean up repo
% -Investigate shoot_template prm
% -Create pars.seg.
% -Modify write to use missing data
% -Add readme
% -Add spm_preproc_img (introduces negative values!)
% -remove gmm.po (test on IXIs)
% -Same results, cluster vs non-cluster
% -Introduce aux
% -Improved registration
% -Labelled resps
% -K.rem,K.keep,K.lkp

addpath(genpath('./code'))
addpath('/cherhome/mbrud/dev/distributed-computing/')

%--------------------------------------------------------------------------
dir_output = '/data-scratch/mbrud/data/build-template/';

%-------------------------------------------------------------------------- 
% K           = 12;
% im          = {};
% % im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-CHROMIS-noneck/',...
% %                Inf,'CT',0,3,'mean',''};
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-healthy-noneck/',...
%                20,'CT',[4 11],3,'mean',''};            
% % im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/IXI-noneck/',...
% %                15,'MRI',3,3,'mean',''}; % 60
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-big-lesions/',...
%                20,'CT',[],3,'mean',''}; % 79 
           
% Uncomment below for testing     
S           = 1;
K           = 6;           
im          = {};
im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/IXI-noneck/',...
               S,'MRI',[],4,'mean',''}; 
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-CHROMIS-noneck/',...
%                S,'CT',[],4,'random',''};

%--------------------------------------------------------------------------
holly = struct;

holly.server.ip      = 'holly';
holly.server.login   = 'mbrud';
holly.server.folder  = ['/scratch/' dir_output(15:end)];
holly.client.folder  = dir_output;

% Uncomment the below to run a regular for-loop
holly.server.ip      = '';
holly.client.workers = 0;

holly.matlab.bin     = '/share/apps/matlab';
holly.matlab.add     = '/home/mbrud/dev/build-template/';

holly.translate      = {'/data-scratch/mbrud/' '/scratch/mbrud/'};
holly.restrict       = 'char';
holly.clean          = true;
holly.verbose        = false;

holly.job.est_mem    = true;
holly.job.batch      = true;
holly.job.mem        = '4G';
holly.job.use_dummy  = true;

%--------------------------------------------------------------------------
pars.print_ll  = true;
pars.print_seg = true;

pars.ml           = false;
pars.do_bf        = true;
pars.do_def       = true;
pars.do_wp        = true;
pars.do_write_res = false;
pars.kmeans_dist  = 'cityblock';

pars.niter = 30;

pars.dir_output = dir_output;

do_show_seg     = true;
do_shrink       = true;
do_show_results = 3;
tol             = 1e-4;

pars.missing_data   = false;
do_avg_template_dim = true;

%--------------------------------------------------------------------------
M               = numel(im);
V               = cell(1,M);
for m=1:M, V{m} = get_V(im{m}); end

%--------------------------------------------------------------------------
pth_template = '';

% Uncomment below to use predefined templates
% pth_template = fullfile(spm('dir'),'tpm','TPM.nii');  
% pars.lkp     = [1 1 2 2 3 3 4 4 5 5 5 6 6];
% pth_template = fullfile(get_pth_dropbox,'/PhD/Data/CB-TPM/BlaiottaTPM.nii');
% pars.lkp     = [1 1 2 2 3 3 4 4 5 5 6 6 7 7];

[pth_template,uniform] = init_template(pth_template,V,K,pars.dir_output,do_avg_template_dim); 

%------------------------------------------------------
[obj,niter,fig,rand_subjs] = init_obj(V,im,pth_template,uniform,do_show_seg,do_show_results,pars);

%-------------------------------------------------------------------------- 
holly = distribute_default(holly);

%--------------------------------------------------------------------------
print_algorithm_progress('started')

L = -Inf;
for iter=1:niter        
               
    % Compute average bias field DC component from pth_obj
    %----------------------------------------------------------------------    
    obj = get_avg_bf_dc(obj);       
    
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
        obj = update_intensity_prior(obj);
    end
       
    if niter>1
        % Some verbose
        %----------------------------------------------------------------------
        if do_show_results>0, plot_ll(fig{5},L); end
        if do_show_results>1, show_template(fig{6},pth_template,obj); end
        if do_show_results>2, show_resp(fig{7},obj,rand_subjs); end      

        tol1 = abs((L(end - 1)*(1 + 10*eps) - L(end))/L(end));    
        fprintf('%i | L = %0.0f | diff = %0.7f | tol = %0.5f\n',iter,L(end),L(end) - L(end - 1),tol1);  
    end
end

print_algorithm_progress('finished')
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
    fprintf('%i | size(otpm) = [%d %d %d] | size(ntpm) = [%d %d %d]\n',iter,od(1),od(2),od(3),nd(1),nd(2),nd(3));
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
function obj = get_avg_bf_dc(obj)
M = numel(obj);
for m=1:M
    S         = numel(obj{m});
    sum_bf_dc = 0;
    cnt_S     = 0;
    for s=1:S                
        if obj{m}{s}.status==0
            sum_bf_dc = sum_bf_dc + obj{m}{s}.bf_dc;
            cnt_S     = cnt_S + 1;
        end
    end
    
    avg_bf_dc = sum_bf_dc/cnt_S;
    
    for s=1:S 
        obj{m}{s}.avg_bf_dc = avg_bf_dc; 
    end
end
%==========================================================================

%==========================================================================
function show_template(fig,pth_template,obj)
M  = numel(obj);
pr = [];
for m=1:M
    if strcmp(obj{m}{1}.modality,'CT')
        pr = obj{m}{1}.gmm.pr;
        
        break
    end
end

set(0,'CurrentFigure',fig);     

Nii = nifti(pth_template);
b   = Nii.dat(:,:,:,:);
d   = size(b);
K   = d(4);                                  
for k=1:K         
    subplot(3,K,k);
    slice = b(:,:,floor(d(3)/2) + 1,k);
    imagesc(slice'); axis image xy off; title(['k=' num2str(k)]); colormap(pink);               
   
    subplot(3,K,K + k);
    slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
    if ~isempty(pr)
        imagesc(slice); axis image xy off; title(['m=' num2str(round(pr.m(1,k),2))]); colormap(pink);   
    else
        imagesc(slice); axis image xy off; colormap(pink);   
    end

    subplot(3,K,2*K + k);
    slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
    imagesc(slice'); axis image xy off; colormap(pink);   
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
            imagesc(img'); axis image xy off; colormap(pink);
            title(['q_{' num2str(m), ',' num2str(s) ',' num2str(k) '}']);
            cnt = cnt + 1;
        end 
    end
end                                  
drawnow
%==========================================================================

%==========================================================================
function print_algorithm_progress(status)
date    = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=================================================\n')  
fprintf('| %s | Algorithm %s\n',date,status)
fprintf('=================================================\n\n')   
%==========================================================================

%==========================================================================
function plot_ll(fig,L)
set(0,'CurrentFigure',fig);                
plot(0:numel(L(3:end)) - 1,L(3:end),'b-','LineWidth',1);   hold on            
plot(0:numel(L(3:end)) - 1,L(3:end),'b.','markersize',10); hold off  
title('ll')
%==========================================================================