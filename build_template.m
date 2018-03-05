function build_template

%--------------------------------------------------------------------------
%
% Qs
% -Is something wrong with shrinking? (not shrinking after iter=1...)
%
% PRIO
% -Add back missing 
%   --problem with convergence of bf (when lkp>1?), probably due to
%     decreasing ll...
%   --include missing ll for ml
%
% -Clean up repo + add readme
% -Labelled resps
% -Improved registration: diffeo + reparameterise + average correct
%
% TODO
% -do_reslice not working
% -impreproc -> overwrite or not
% -Init parameters in a better way
% -Test including MRA images (requires missing data to work well)
% -DICOM convert
% -Add denoising
% -Add superres (spm_imquality)
% -Introduce aux
% -Use K - 1 classes for template?
% -Investigate shoot_template prm
% -K.rem,K.keep,K.lkp
%
% CLUSTER
% -Reintroduce split?
% -Same results, cluster vs non-cluster
% -If same result, try to sum to one gr and H
%
%--------------------------------------------------------------------------

addpath(genpath('./code'))
addpath('/cherhome/mbrud/dev/distributed-computing')
addpath('/cherhome/mbrud/dev/auxiliary-functions')

TEST_LEVEL = 0;

pars              = [];
pars.name         = 'CT-lesion-holly';
pars.dir_output   = '/data-scratch/mbrud/data/';
pars.dir_template = '/data/mbrud/templates/';
% pars.dir_output   = '/home/mbrud/Data/temp-data/';
% pars.dir_template = '/home/mbrud/Data/template/';

%--------------------------------------------------------------------------
% Set directories to input images as well as where to write output
%--------------------------------------------------------------------------

im = {};       
if TEST_LEVEL
    if TEST_LEVEL<3, S = 1;
    else             S = 8;
    end
    
%     pars.K      = 6;  
%     im{end + 1} = {'/data/mbrud/images/MRI/IXI-T1-preproc-rn/',...
%                    S,'MRI',[],3,'','',''};    

    pars.K      = 8;            
    im{end + 1} = {'/data/mbrud/images/CT/CHROMIS-preproc-rn-ss/',...
                   8,'CT',[],3,'','',''};  
    im{end + 1} = {'/data/mbrud/images/CT/healthy-preproc-rn-ss/',...
                   8,'CT',[7],3,'','',''};      
               
    if TEST_LEVEL<4, pars.name = 'test-local';
    else             pars.name = 'test-holly';   
    end               
else    
    pars.K      = 8;            
    im{end + 1} = {'/data/mbrud/images/CT/CHROMIS-preproc-rn-ss/',...
                   200,'CT',[],3,'','',''};  
    im{end + 1} = {'/data/mbrud/images/CT/healthy-preproc-rn-ss/',...
                   Inf,'CT',[7],3,'','',''};                   
end

pars = append_dir(pars);

%--------------------------------------------------------------------------
% Set-up parallel options
%--------------------------------------------------------------------------

holly               = struct;
holly.server.ip     = 'holly';
holly.server.login  = 'mbrud';
holly.server.folder = fullfile('/scratch',pars.dir_output(15:end),'cluster');
holly.client.folder = fullfile(pars.dir_output,'cluster');
holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = '/home/mbrud/dev/build-template';
holly.matlab.add    = '/home/mbrud/dev/auxiliary-functions';
holly.translate     = {'/data-scratch/mbrud/' '/scratch/mbrud/'};
holly.restrict      = 'char';
holly.clean         = false;
holly.clean_init    = true;
holly.verbose       = false;
holly.job.est_mem   = true;
holly.job.batch     = true;
holly.job.mem       = '6G';
holly.job.use_dummy = true;

if TEST_LEVEL>0 && TEST_LEVEL<4
    holly.server.ip  = '';   
    
    if TEST_LEVEL<3, holly.client.workers = 0;
    else             holly.client.workers = Inf;
    end
end

% holly.server.ip      = '';   
% holly.client.workers = Inf;

holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Set parameters
%--------------------------------------------------------------------------

pars.do_segment          = true;
pars.do_preproc          = false;
pars.niter               = 50;
pars.tol                 = 1e-4;
pars.verbose             = 3;
pars.pth_template        = '';
% pars.pth_template        = '/home/mbrud/Dropbox/PhD/Data/log-template/logTPM.nii';
% pars.pth_template        = fullfile(get_pth_dropbox,'/PhD/Data/CB-TPM/BlaiottaTPM.nii');
pars.vx_tpm              = 1.5;
pars.sparam              = [0.01 2 0]; 
if TEST_LEVEL==2
    pars.do_preproc      = true;
    pars.do_segment      = false;
end

pars.preproc.do_rem_corrupted  = false;
pars.preproc.tol_dist          = 4;
pars.preproc.tol_vx            = 5;
pars.preproc.verbose           = false;
pars.preproc.coreg_and_reslice = true;
pars.preproc.do_reslice        = true;
pars.preproc.realign2mni       = true;
pars.preproc.crop              = true;
pars.preproc.rem_neck          = true;
pars.preproc.skull_strip       = true;

pars.segment.do_ml           = false;
pars.segment.do_bf           = true;
pars.segment.do_def          = true;
pars.segment.do_wp           = true;
pars.segment.do_write_res    = false;
pars.segment.do_push_resp    = false;
pars.segment.do_missing_data = false;
pars.segment.do_old_segment  = false;
pars.segment.kmeans_dist     = 'cityblock';
pars.segment.print_ll        = false;
pars.segment.print_seg       = false;    
pars.segment.verbose         = false;   
pars.segment.trunc_ct        = [-Inf Inf];   
pars.segment.nlkp            = 2;
pars.segment.lkp             = reshape(repmat(1:pars.K,2,1),1,[]);

if TEST_LEVEL==1
    pars.segment.verbose   = true;
    pars.segment.print_ll  = true;
    pars.segment.print_seg = true;
end

%--------------------------------------------------------------------------
% Init
%--------------------------------------------------------------------------

V = cell(1,numel(im));
for m=1:numel(im) 
    V{m} = read_images(im{m},pars); 
end

% for m=1:M, 
%     V{m} = read_labels(V{m},im{m}); 
% end

pars = init_template(V,pars); 

[obj,pars,fig,rand_subjs] = init_obj(V,im,pars);

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
        obj = modify_obj(obj,pars,iter);
    end       
    
    % Segment a bunch of subjects 
    %----------------------------------------------------------------------
    [obj,ix]    = unfold_cell(obj,2);
    [holly,obj] = distribute(holly,'process_subject','inplace',obj,fig);
    obj         = fold_cell(obj,ix);

    % Check if any subjects have status~=0
    %----------------------------------------------------------------------
    print_jobs_failed(obj);
    
    if pars.niter>1
        % Update template
        %------------------------------------------------------------------
        L = update_template(L,obj,pars.sparam,iter);                              
    end        
    
    if pars.niter>1
        % Automatically decrease the size of the template based on bounding
        % boxes computed for each pushed responsibility image
        %------------------------------------------------------------------
        shrink_template(obj,iter);    
    end
    
    if pars.niter>1
        % Update Gaussian-Wishart hyper-parameters
        %------------------------------------------------------------------
        obj = update_intensity_prior(obj,iter);
    end
       
    % Save object file
    %----------------------------------------------------------------------
    save(fullfile(pars.dir_template,'obj.mat'),'obj')
    
    if pars.niter>1
        % Some verbose
        %----------------------------------------------------------------------
        if pars.verbose>0, plot_ll(fig{5},L); end
        if pars.verbose>1, show_template(fig{6},pars.pth_template); end
        if pars.verbose>2, show_resp(fig{7},obj,rand_subjs); end      

        d = abs((L(end - 1)*(1 + 10*eps) - L(end))/L(end));    
        fprintf('%2d | L = %0.0f | d = %0.5f\n',iter,L(end),d);  
        
        if d<pars.tol
           break 
        end
    end
end

print_algorithm_progress('finished',iter);
%==========================================================================

%==========================================================================
function obj = modify_obj(obj,pars,iter)
K    = pars.K;
nlkp = pars.segment.nlkp;

M = numel(obj);    
for m=1:M
    S         = numel(obj{m});    
    sum_bf_dc = 0;
    cnt_S     = 0;
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==1 
            obj{m}{s}.lkp          = 1:K;
            obj{m}{s}.do_def       = false;
            obj{m}{s}.do_bf        = false;
            obj{m}{s}.do_push_resp = true;
            obj{m}{s}.do_write_res = false;
            
            obj{m}{s}.niter  = 1;                 
            obj{m}{s}.nsubit = 1;
            obj{m}{s}.nitgmm = 1;
        end

        if iter==2
            obj{m}{s}.lkp     = reshape(repmat(1:K,nlkp,1),1,[]);
            obj{m}{s}.nsubit  = 8;
            obj{m}{s}.nitgmm  = 20;  
            obj{m}{s}.do_bf   = true;                
            obj{m}{s}.uniform = false;            
            obj{m}{s}.do_def  = true; 
        end

        if iter>=2                  
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
function pars = append_dir(pars)
pth               = fileparts(pars.dir_output);
pars.dir_output   = fullfile(pth,['build-template-' pars.name]);
pth               = fileparts(pars.dir_template);
pars.dir_template = fullfile(pth,pars.name);
%==========================================================================