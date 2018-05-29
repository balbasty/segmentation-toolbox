function build_template(pars,test_level)
% Builds a probabilistic template from populations of brain images
% FORMAT build_template(pars,test_level)
% pars - Parameter struct (for more info see the function pars_default)
% test_level - For testing and debugging (0=no testing, 1=for, 2=parfor, 3=holly)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<1 
    pars = '/data/mbrud/jobs/segmentation-toolbox/ROB-IXI_3d.json'; 
end
if nargin<2, test_level = 3; end

spm_field('boundary',1); % BOUND_CIRCULANT 0, BOUND_NEUMANN 1

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
pth_distributed_toolbox = '/data/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/data/mbrud/dev/auxiliary-functions';

holly_server_login   = 'mbrud';
holly_matlab_add_src = '/data/mbrud/dev/segmentation-toolbox';
holly_matlab_add_aux = pth_auxiliary_functions;

% addpath
%--------------------------------------------------------------------------

addpath(genpath('./code'))
addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

%--------------------------------------------------------------------------
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
holly.server.folder = ['/' holly.client.folder(7:end)];
holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = holly_matlab_add_src;
holly.matlab.add    = holly_matlab_add_aux;
holly.translate     = {'/data-scratch/','/scratch/'};
holly.restrict      = 'char';
holly.clean         = false;
holly.clean_init    = true;
holly.verbose       = true;
holly.job.batch     = true;
holly.job.mem       = '6G';
holly.job.est_mem   = false;
holly.job.use_dummy = false;
holly.mode          = 'qsub';

if     test_level==1, holly.mode = 'for';
elseif test_level==2, holly.mode = 'parfor';
end

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
    
    if pars.do_template
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
    
    if pars.do_template
        % Update template
        %------------------------------------------------------------------
        L = update_template(L,obj,pars,iter);                              
    end        
                
    if pars.do_template
        % Automatically decrease the size of the template based on bounding
        % boxes computed for each pushed responsibility image
        %------------------------------------------------------------------
        shrink_template(obj,iter);    
    end
    
    if pars.do_template && pars.crop_template==iter
        % Crop template to size of the default SPM template
        %------------------------------------------------------------------
        crop_template(pars.pth_template,iter);
    end
    
    if pars.do_template && pars.mrf
        % Use a MRF cleanup procedure
        %------------------------------------------------------------------
        mrf_template(pars.pth_template);
    end
    
    %----------------------------------------------------------------------
    % Intensity specific
    % Here is where the hyper-parameters of the VB-GMM are updated
    %----------------------------------------------------------------------
    
    if pars.do_template
        % For mean correcting bias field (also updates posteriors)
        %------------------------------------------------------------------
        obj = modify_bf_dc(obj,iter);  
    end  
    
    if pars.do_template
        % Update Gaussian-Wishart hyper-parameters
        %------------------------------------------------------------------
        obj = update_intensity_hp2(obj);
    end
           
    if pars.do_template
        % Show some verbose        
        show_results(pars,obj,L,iter);
                      
        % Save obj struct
        save(fullfile(pars.dir_template,'obj.mat'),'obj');
       
        d = abs((L(end - 1)*(1 + 10*eps) - L(end))/L(end));    
        fprintf('%2d | L = %0.0f | d = %0.5f\n',iter,L(end),d);  
        
%         if d<pars.tol && iter>10
%            break 
%         end
    end
end

print_algorithm_progress('finished',iter);
%==========================================================================

%==========================================================================
function show_results(pars,obj,L,iter)
fig            = pars.fig;
rand_subjs     = pars.rand_subjs;
pth_template   = pars.pth_template;
verbose        = pars.verbose;
dir_animations = pars.dir_animations;

% Model specific verbose
%--------------------------------------------------------------------------
if verbose>0, 
    show_L(fig{5},L,iter); 
end
if verbose>1, 
    show_template(fig{6},pth_template,dir_animations,iter); 
    fname = fullfile(dir_animations,'template.gif');
    write2gif(fig{6},fname,iter);
end

% Subject specific verbose
%--------------------------------------------------------------------------
if verbose>2, 
    show_resp(fig{7},obj,rand_subjs,iter); 
    fname = fullfile(dir_animations,'responsibilities.gif');
    write2gif(fig{7},fname,iter);
end      
if verbose>2, 
    show_intensity_hp(fig{8},obj); 
end 
if verbose>2, 
    show_subj_ll(fig{9},obj,rand_subjs); 
end
if verbose>2, 
    show_def(fig{10},obj,rand_subjs); 
    fname = fullfile(dir_animations,'deformations.gif');
    write2gif(fig{10},fname,iter);
end                 
if verbose>3, 
    show_bf(fig{11},obj,rand_subjs); 
    fname = fullfile(dir_animations,'bias-fields.gif');
    write2gif(fig{11},fname,iter);
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
function show_L(fig,L,iter)
set(0,'CurrentFigure',fig);                
plot(0:numel(L(3:end)) - 1,L(3:end),'b-','LineWidth',1);   hold on            
plot(0:numel(L(3:end)) - 1,L(3:end),'b.','markersize',10); hold off  
title(['L (iter=' num2str(iter) ')']) 
set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]); 
%==========================================================================