function build_template(pars,test_level)
if nargin<1, pars       = struct; end
if nargin<2, test_level = 2; end % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 8 subjects (holly)
if ~isfield(pars,'dir_output')
%     pars.dir_output = '/data/mbrud/tmp-build-tpm/';
    pars.dir_output = '/home/mbrud/Data/temp-segmentation-toolbox';    
end
if ~isfield(pars,'name')
    pars.name = 'segmentation-toolbox';  
end

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------

m = 0;

% Basic test
% -----------------
m = m + 1;

pars.dat{m}.dir_data = '/home/mbrud/Dropbox/PhD/Data/IXI-test/2d_IXI-T1T2PD_preproc-ra-cr-rn-reg-res-vx';
% pars.dat{m}.dir_data = '/home/mbrud/Dropbox/PhD/Data/IXI-test/IXI-T1T2PD_preproc-ra-cr-rn-reg-res-vx';
pars.dat{m}.S = 8;
pars.dat{m}.segment.samp = 1;

% CT
% -----------------
% pars.K = 20;
% 
% m = m + 1;
% 
% pars.dat{m}.dir_data = '/data/mbrud/images/CT/test-labels_preproc-ra-cr-rn-vx';
% % pars.dat{m}.dir_data = '/data/mbrud/images/CT/2d_test-labels_preproc-ra-cr-rn-vx-1';
% pars.dat{m}.modality = 'CT';
% pars.dat{m}.S = Inf;
% pars.dat{m}.segment.samp = 1;
% 
% pars.dat{m}.segment.kmeans_hist = true;
% pars.dat{m}.segment.do_bf = false;
% 
% pars.dat{m}.segment.lkp.lab = zeros(1,pars.K);
% pars.dat{m}.segment.lkp.lab(15) = 1;
% pars.dat{m}.segment.wp_lab = 0.5;

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
        mrf_template(pars.pth_template);
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
        %----------------------------------
% FORMAT nii2subject(dir_in,dir_out)--------------------------------
        obj = update_intensity_prior(obj,iter);
    end
       
    % Save obj structs
    %------------------------------------------------------------------
    save(fullfile(pars.dir_template,'obj.mat'),'obj');
    
    if pars.niter>1
        % Some verbose
        %------------------------------------------------------------------
        if pars.verbose>0, plot_ll(pars.fig{5},L); end
        if pars.verbose>1, show_template(pars.fig{6},pars.pth_template); end
        if pars.verbose>2, show_resp(pars.fig{7},obj,pars.rand_subjs); end      
        if pars.verbose>3, show_def(pars.fig{8},obj,pars.rand_subjs); end 

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