function build_template

%--------------------------------------------------------------------------
%
% Qs
% -Is something wrong with shrinking? (not shrinking after iter=1...)
%
% PRIO
% -labels
%   --FSL_SEG -> .img
%   --PROCESSED -> OAS1_0001_MR1_mpr_n4_anon_111_t88_gfc
%   --add to preprocess as well
% -Add back missing 
%   --problem with convergence of bf (when lkp>1?), probably due to
%     decreasing ll...
%   --include missing ll for ml
% -Clean up repo + add readme
% -Improved registration: diffeo + reparameterise + average correct
%
% TODO
% -do_reslice not working
% -impreproc -> overwrite or not
% -Test including MRA images (requires missing data to work well)
% -DICOM convert
% -Add denoising
% -Add superres (spm_imquality)
% -Introduce aux
% -Use K - 1 classes for template?
% -Investigate shoot_template prm
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

test_level = 0;

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------

pars.name         = 'CHROMIS';
pars.dir_output   = '/data/mbrud/data-seg';
% pars.dir_output   = '/data-scratch/mbrud/data-seg/';
% pars.dir_output   = '/home/mbrud/Data/data-seg';
pars.dat          = {};

% m = 1;
% pars.dat{m}.dir_data = '/data/mbrud/images/MRI/IXI-T1T2PD/';
% pars.dat{m}.do_preproc = true;
% pars.dat{m}.do_rem_corrupted = false;
% pars.dat{m}.S = Inf;
% pars.do_segment = false;

% m = 1;
% pars.dat{m}.dir_data = '/home/mbrud/Data/MR-and-labels/';
% pars.dat{m}.segment.use_labels = true;
% % pars.dat{m}.segment.wp_lab     = eps;
% pars.dat{m}.segment.do_push_resp = true;
% pars.dat{m}.segment.ml = true;
% pars.dat{m}.segment.lkp.lab = [3 2 1 0 0 0];
% % % pars.dat{m}.segment.lkp.rem = [6];

pars.K = 9;  
m = 1;
pars.dat{m}.dir_data           = '/data/mbrud/images/CT/CHROMIS-preproc-rn-ss/';
pars.dat{m}.S                  = Inf;
pars.dat{m}.modality           = 'CT';
pars.dat{m}.segment.samp       = 2;
pars.dat{m}.segment.init_clust = 'total';

m = 2;
pars.dat{m}.dir_data           = '/data/mbrud/images/CT/healthy-preproc-rn-ss/';
pars.dat{m}.S                  = Inf;
pars.dat{m}.modality           = 'CT';
pars.dat{m}.segment.samp       = 2;
pars.dat{m}.segment.init_clust = 'total';
pars.dat{m}.segment.lkp.rem    = [4];        

% pars.pth_template        = fullfile(get_pth_dropbox,'PhD/Data/log-template/logTPM.nii');
% pars.pth_template        = fullfile(get_pth_dropbox,'PhD/Data/log-template/logBlaiottaTPM.nii');

pars = pars_default(pars,test_level);

%--------------------------------------------------------------------------
% Set distribute package parameters
%--------------------------------------------------------------------------

holly               = struct;
holly.server.ip     = 'holly';
holly.server.login  = 'mbrud';
holly.client.folder = fullfile(pars.dir_output,'cluster');
holly.server.folder = holly.client.folder;
% holly.server.folder = fullfile('/scratch',holly.client.folder(15:end));
holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = '/home/mbrud/dev/build-template';
holly.matlab.add    = '/home/mbrud/dev/auxiliary-functions';
% holly.translate     = {'/data-scratch/mbrud/' '/scratch/mbrud/'};
holly.restrict      = 'char';
holly.clean         = false;
holly.clean_init    = true;
holly.verbose       = false;
holly.job.mem       = '6G';
holly.job.use_dummy = true;

if     test_level==1, holly.server.ip  = ''; holly.client.workers = 0;
elseif test_level==2  holly.server.ip  = ''; holly.client.workers = Inf;
end
% holly.server.ip      = ''; holly.client.workers = Inf;
% holly.server.ip      = ''; holly.client.workers = 0;

holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Initialise algorithm
%--------------------------------------------------------------------------

pars       = read_images(pars); 
pars       = init_template(pars); 
[obj,pars] = init_obj(pars);

%--------------------------------------------------------------------------
% Start the algorithm
%--------------------------------------------------------------------------

print_algorithm_progress('started');

L = -Inf;
for iter=1:pars.niter        
    
    if pars.niter>1 && pars.sum_temp_der
        clear_temp_der(obj);
    end
    
    if pars.niter>1
        % Some parameters of the obj struct are changed depending on iteration 
        % number (only for building templates)
        %----------------------------------------------------------------------    
        obj = modify_obj(obj,iter);
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
        if pars.verbose>0, plot_ll(pars.fig{5},L); end
        if pars.verbose>1, show_template(pars.fig{6},pars); end
        if pars.verbose>2, show_resp(pars.fig{7},obj,pars); end      

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
            obj{m}{s}.do_bf   = true;                
            obj{m}{s}.uniform = false;                        
        end

        if iter>=3            
            obj{m}{s}.reg    = obj{m}{s}.reg0;            
            scal             = 2^max(12 - iter,0);                  
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
function show_template(fig,pars)
set(0,'CurrentFigure',fig);     

pth_template = pars.pth_template;

Nii = nifti(pth_template);
b   = exp(Nii.dat(:,:,:,:));
b   = bsxfun(@rdivide,b,sum(b,4));
d   = size(b);
K   = d(4);                                  
for k=1:K         
    subplot(3,K,k);
    slice = b(:,:,floor(d(3)/2) + 1,k);
    imagesc(slice'); axis off xy; title(['k=' num2str(k)]); colormap(gray);               
   
    subplot(3,K,K + k);
    slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
    imagesc(slice); axis off xy; colormap(gray);   

    subplot(3,K,2*K + k);
    slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
    imagesc(slice'); axis off xy; colormap(gray);   
end 
drawnow
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
            zix = floor(dm(3)/2) + 1;
            img = img(:,:,zix);
        
            subplot(M*numel(rand_subjs{1}),K,cnt);
            imagesc(img'); axis off xy; colormap(gray);
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
function clear_temp_der(obj)
V    = spm_vol(obj{1}{1}.pth_template); 
d1   = [V(1).dim numel(V)];
d_gr = [d1(1:3),d1(4) - 1];
d_H  = [d1(1:3) round(((d1(4) - 1)*d1(4))/2)];

create_nii(obj{1}{1}.pth_sgr,zeros(d_gr,'single'),eye(4),obj{1}{1}.dt,'sgr'); 
create_nii(obj{1}{1}.pth_sH,zeros(d_H,'single'),eye(4),obj{1}{1}.dt,'sH'); 
%==========================================================================