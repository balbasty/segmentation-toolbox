function build_template(pars,test_level)
% Builds a probabilistic template from populations of brain images
% FORMAT build_template(pars,test_level)
% pars - Parameter struct (for more info see the function pars_default)
% test_level - For testing and debugging (0=no testing, 1=for, 2=parfor, 3=holly)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, test_level = 0; end

if test_level==0 || test_level==3
    dir_pars = '/data/mbrud/pars/segmentation-toolbox/'; 
else
    dir_pars = '/home/mbrud/Dropbox/PhD/Data/pars/segmentation-toolbox/'; 
end

if nargin<1 
    pars = [dir_pars 'CROMIS-and-healthy-3d-vx.json']; 
end

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
if test_level==0 || test_level==3
    pth_distributed_toolbox = '/cherhome/mbrud/dev/distributed-computing';
    pth_auxiliary_functions = '/cherhome/mbrud/dev/auxiliary-functions';
else
    pth_distributed_toolbox = '/home/mbrud/dev/WTCN-computational-anatomy-group/distributed-computing';
    pth_auxiliary_functions = '/home/mbrud/dev/WTCN-computational-anatomy-group/auxiliary-functions';
end

holly_server_login   = 'mbrud';
holly_matlab_add_src = '/home/mbrud/dev/segmentation-toolbox';
holly_matlab_add_aux = '/home/mbrud/dev/auxiliary-functions';

% addpath
%--------------------------------------------------------------------------

addpath(genpath('./code'))
addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

%--------------------------------------------------------------------------\
% Set algorithm parameters
%--------------------------------------------------------------------------
pars = segment_default(pars,test_level); 

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
holly.verbose       = true;
holly.job.batch     = true;
holly.job.mem       = '4G';
holly.job.est_mem   = true;
holly.job.use_dummy = false;

if     test_level==1, holly.server.ip  = ''; holly.client.workers = 0;
elseif test_level==2, holly.server.ip  = ''; holly.client.workers = Inf;
end
% holly.server.ip = ''; holly.client.workers = Inf;
% holly.server.ip = ''; holly.client.workers = 0;

holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Initialise algorithm
%--------------------------------------------------------------------------
pars = read_images_segment(pars); 
pars = init_template(pars); 

% Only if CT data, initialise the GMM parameters by fitting a GMM to an
% accumulated histogram of image intensities.
pars = init_ct_gmm(pars);

[obj,pars] = init_obj_segment(pars);

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
    
    %----------------------------------------------------------------------
    % Segmentation part
    % Here is where the VB-GMM, bias field, and registration parameters are
    % updated
    %----------------------------------------------------------------------
    
    [obj,ix]    = unfold_cell(obj,2);
    [holly,obj] = distribute(holly,'process_subject_segment','inplace',obj,pars.fig);
    obj         = fold_cell(obj,ix);
    
    % Check if any subjects have status~=0
    print_jobs_failed(obj);
           
    %----------------------------------------------------------------------
    % Template specific
    % Here is where the log-template is updated
    %----------------------------------------------------------------------
    
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
    
    %----------------------------------------------------------------------
    % Intensity specific
    % Here is where the hyper-parameters of the VB-GMM are updated
    %----------------------------------------------------------------------
    
    if pars.niter>1
        % For mean correcting bias field (also updated posteriors)
        %------------------------------------------------------------------
        obj = calc_avg_dc(obj);  
    end  
    
    if pars.niter>1
        % Update Gaussian-Wishart hyper-parameters
        %------------------------------------------------------------------
        obj = update_intensity_hp2(obj,iter,pars);
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
        if pars.verbose>2, show_intensity_hp(pars.fig{8},obj); end 
        if pars.verbose>2, show_subj_ll(pars.fig{9},obj,pars.rand_subjs); end
        if pars.verbose>2, show_def(pars.fig{10},obj,pars.rand_subjs); end                 
        if pars.verbose>3, show_bf(pars.fig{11},obj,pars.rand_subjs); end 

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