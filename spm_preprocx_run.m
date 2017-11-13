function spm_preprocx_run

close all;

addpath(genpath('./code'))

% -148243

%--------------------------------------------------------------------------
% Define data cell array, which should contain the following:
% {'pth_to_images',num_subjects,'CT'/'MRI','healthy'/'non-healthy','pth_to_labels'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.

% %-----------------
% % Russel Sq. House
% S        = Inf;
% Kb       = 16;
% imobj{1} = {'./tmp/',S,'CT','healthy',''};

%-----------------
% % Aged 2D MRI+CT
% S        = 256; % Number of subjects
% Kb       = 16; % Number of classes (if a template is used, then Kb will be set to the number of classes in that template)
% imobj{1} = {'/home/mbrud/Data/2D-Data/OASIS-Longitudinal-2D/',S,'MRI','healthy',''};
% imobj{2} = {'/home/mbrud/Data/2D-Data/CT-aged-2D/',S,'CT','healthy',''};

%-----------------
% IXI 2D
S        = 1;
Kb       = 6;
imobj{1} = {'/home/brudfors/Dropbox/PhD/Data/2D-Data/IXI-2D/',S,'MRI','healthy',''};

% %-----------------
% % CT 2D
% S        = Inf; % Number of subjects
% Kb       = 16; % Number of classes (if a template is used, then Kb will be set to the number of classes in that template)
% imobj{1} = {'/home/brudfors/Dropbox/PhD/Data/2D-Data/CT-hemorrhage-2D/',S,'CT','healthy',''};
% % imobj{1} = {'/home/mbrud/Data/2D-Data/CT-aged-2D/',S,'CT','healthy',''};

%--------------------------------------------------------------------------
% Preprocessing options
preproc.is_DICOM      = 0; % If input images are DICOM, converts DICOM to Nifti % TODO (work in progress)
preproc.rem_corrupted = 1; % Try to remove CT images that are corrupted (e.g. bone windowed)
preproc.do_preproc    = 0; % Do preprocessing on input images
preproc.realign       = 1; % Realign to MNI space
preproc.crop          = 1; % Remove data outside of head
preproc.crop_neck     = 0; % Remove neck (the spine, etc.)
preproc.denoise       = 1; % Denoise (only done if data is CT)

%--------------------------------------------------------------------------
% Run the algorithm in parallel
runpar = Inf;

%--------------------------------------------------------------------------
% The distance (mm) between samples (for sub-sampling input data--improves speed)
samp = 2;

%--------------------------------------------------------------------------
% Segmentation parameters
nlkp      = 1; % Number of gaussians per tissue
use_vbmog = 0; % Use a variational Bayesian mixture model

%--------------------------------------------------------------------------
% What estimates to perform
domiss = 0; % Handle missing data
dobias = 1; % Bias field
doaff  = 1; % Affine registration
dodef  = 1; % Non-linear registration
dopr   = 1; % Intensity priors
dotpm  = 1; % Template

%--------------------------------------------------------------------------
% Different number of iterations and stopping tolerance of algorithm
nitermain = 200;
tolmain   = 1e-4;

niter     = 1;
niter1    = 20;
nsubitmog = 20;
nsubitbf  = 1;
nitdef    = 1;

%--------------------------------------------------------------------------
% Regularisation for estimating deformations
rparam = [0 0.001 0.5 0.05 0.2]*0.1;

%--------------------------------------------------------------------------
% Template options
pth_logTPM = '';   % Path to existing template (set to '' for estimating a template, or get_spm_TPM for using the default SPM one)
vx_TPM     = 1.5;  % Voxel size of template to be estimated
deg        = 2;    % Degree of interpolation when sampling template
tiny       = 1e-4; % Strength of Dirichlet prior used in template construction
fwhm_TPM   = 0.1;  % Ad hoc smoothing of template
mrf        = 0;    % Use a MRF cleanup procedure

%--------------------------------------------------------------------------
% For debugging
verbose = 1;
figix   = 1;

%--------------------------------------------------------------------------
% Make some folders
dir_data = './data'; % for all algorithm data
if ~exist(dir_data,'dir'), mkdir(dir_data); end
dir_res = fullfile(dir_data,'results'); % for algorithm results
if exist(dir_res,'dir'), rmdir(dir_res,'s'); end; mkdir(dir_res);

%--------------------------------------------------------------------------
if ~isempty(pth_logTPM) || S==1   
    runpar    = 0;
    dopr      = 0;
    dotpm     = 0;
    nitermain = 1;
    niter     = 30;
    niter1    = 8;
    nsubitmog = 20;
    nsubitbf  = 1;
    nitdef    = 3;
end

%==========================================================================
%% Load and process image data
[V,M,S,N,labels] = load_and_process_images(imobj,preproc,dir_data,runpar);

%==========================================================================
%% Initialise template
if isempty(pth_logTPM)
    pth_logTPM = init_logTPM(V,Kb,vx_TPM,dir_res);
    use_tpm    = false;    
    uniform    = true;
else    
    Nii     = nifti(pth_logTPM);
    Kb      = size(Nii.dat(:,:,:,:),4);
    use_tpm = true;    
    uniform = false;    
    clear Nii
end

%==========================================================================
%% Initialise algorithm i/o   
fig = cell(4,1);
if verbose    
    if ~runpar
        spm_figure('Create','Interactive');
        figure(figix)
    
        for i=1:size(fig,1)              
            fig{i} = figure(figix + i);
            clf(fig{i})
        end 
    end

    fig_TPM = figure(figix + 5);
    fig_L   = figure(figix + 6);
    
    try
        distFig; 
    catch
        warning('distFig not available')
    end
    drawnow;
end  
    
dir_Twarp = fullfile(dir_data,'Twarp');
if exist(dir_Twarp,'dir'), rmdir(dir_Twarp,'s'); end; mkdir(dir_Twarp);

obj = cell(1,M);
for m=1:M    
    obj{m} = cell(1,S(m));
    for s=1:S(m)
        obj{m}{s}.image    = V{m}{s};
        obj{m}{s}.biasfwhm = 60*ones(1,N(m));
        obj{m}{s}.biasreg  = 1e-3*ones(1,N(m));       

        obj{m}{s}.use_vbmog = use_vbmog;
        obj{m}{s}.use_tpm   = use_tpm;

        obj{m}{s}.lkp  = repelem(1:Kb,nlkp);
        obj{m}{s}.nlkp = nlkp;

        obj{m}{s}.Affine  = eye(4);
        obj{m}{s}.reg     = rparam;
        obj{m}{s}.samp    = samp;
        obj{m}{s}.fwhm    = 0;
        obj{m}{s}.verbose = verbose;

        obj{m}{s}.domiss         = domiss;
        obj{m}{s}.dobias         = dobias;
        obj{m}{s}.dodef0         = dodef; 
        obj{m}{s}.dotpm          = dotpm;         
        obj{m}{s}.doaff          = doaff;
        obj{m}{s}.clear_mog_pars = 1;
        
        obj{m}{s}.niter      = niter;
        obj{m}{s}.niter1     = niter1;
        obj{m}{s}.nsubitmog  = nsubitmog;
        obj{m}{s}.nsubitbf   = nsubitbf;
        obj{m}{s}.nitdef     = nitdef;

        obj{m}{s}.descrip = imobj{m}{3};
        obj{m}{s}.healthy = imobj{m}{4}; 
        obj{m}{s}.labels  = labels{m}{s}; 
        
        obj{m}{s}.fig = fig; 
              
        if S(m)>1
            % Allocate initial deformation field (only for multiple subjects) 
            [obj{m}{s}.pth_Twarp] = alloc_Twarp(obj{m}{s},m,s,dir_Twarp);
        end
    end
end

%==========================================================================
%% Run algorithm
fprintf('==============================================\n')   
fprintf('Algorithm started (')
fprintf(datestr(now,'mmmm dd, yyyy HH:MM:SS'))
fprintf(')\n')
fprintf('==============================================\n')   

Nm = 0;    % Total number of voxels in all images
L  = -Inf; % Lower bound of complete model
for iter=1:nitermain
    fprintf('iter=%d======================\n',iter);  

    % Update the template--------------------------------------------------                     
    if dotpm && iter>1   
        uniform = false;
        
        Nii = nifti(pth_logTPM);
        dm  = size(Nii.dat(:,:,:));
        if numel(dm)==2
            dm(3) = 1;
        end
        
        % Smooth template
        [munum,muden] = smooth_template(munum,muden,dm,fwhm_TPM);

        % Update template
        logmu            = log(munum./muden + tiny);
        logmu            = reshape(logmu',[dm Kb]);                                          
        Nii.dat(:,:,:,:) = logmu;        
        clear munum muden
        
        % Use a MRF cleanup procedure
        if mrf
            logtpm = spm_load_logpriors8(pth_logTPM,tiny,deg,uniform);
            
            phi = double(identity(logtpm.d));
            b   = spm_sample_logpriors8(logtpm,phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3));
            clear phi
        
            Q = zeros([logtpm.d(1:3),Kb],'single');
            for k=1:Kb
               Q(:,:,:,k) = b{k};
            end
            clear b
            
            P        = zeros([logtpm.d(1:3),Kb],'uint8');
            nmrf_its = 1;           
            T        = 1;
            G        = ones([Kb,1],'single')*T;
            vx2      = 1./single(sqrt(sum(logtpm.M(1:3,1:3).^2)));
            for i=1:nmrf_its
                spm_mrf(P,Q,G,vx2);
            end             
            clear Q logtpm
            
            P                = double(P)/255;
            Nii.dat(:,:,:,:) = log(P + eps*eps);
            clear P
        end   
        clear Nii logmu
    end                 
    
    % Load template--------------------------------------------------------
    logtpm = spm_load_logpriors8(pth_logTPM,tiny,deg,uniform);
    
    % Visualise template---------------------------------------------------
    if dotpm && iter>1 && verbose        
        phi = double(identity(logtpm.d));
        b   = spm_sample_logpriors8(logtpm,phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3));
        clear phi

        K1 = floor(sqrt(Kb)); K2 = ceil(Kb/K1); 
        set(0,'CurrentFigure',fig_TPM);                                        
        for i=1:Kb    
            subplot(K1,K2,i);
            imagesc(b{i}(:,:,floor(logtpm.d(3)/2) + 1)'); axis image xy off; colormap(gray);
        end 
        clear b
        drawnow
    end     
    
    % Update subject specific parameters-----------------------------------
    munum = single(0); 
    muden = single(0);  
    ll    = 0;
    for m=1:M
        obj_m = obj{m};
        for s=1:S
%         parfor (s=1:S(m),runpar)    
            fprintf('iter=%d, m=%d, s=%d\n',iter,m,s);  
            
            [obj_m{s},munum1,muden1] = update_subject_pars(obj_m{s},logtpm,iter);
            munum                    = munum + munum1;
            muden                    = muden + muden1;            
            ll                       = ll + obj_m{s}.ll;
            if iter==1
               Nm                    = Nm + obj_m{s}.nm;
            end            
        end
        obj{m} = obj_m;
        clear obj_m
    end
    L = [L,ll];

    % Plot lower bound-----------------------------------------------------
    if verbose        
        set(0,'CurrentFigure',fig_L);                
        plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
        plot(0:numel(L) - 1,L,'b.','markersize',10); hold off                        
    end
       
    % Update intensity prior-----------------------------------------------
    if dopr && use_vbmog && iter>1
        for m=1:M
            pr = update_intensity_prior(obj{m});
            for s=1:S(m)
                obj{m}{s}.pr = pr;
            end      
            save(fullfile(dir_res,['pr_m' num2str(m) '.mat']),'pr'); 
        end
              
        % Update subject specific parameters-----------------------------------        
        ll = 0;
        for m=1:M
            obj_m = obj{m};
    %         for s=1:S
            parfor (s=1:S(m),runpar)    
                obj_m{s} = update_subject_pars(obj_m{s},logtpm,iter);
                ll       = ll + obj_m{s}.ll; % Sum up likelihoods over all subjects                        
            end
            obj{m} = obj_m;
            clear obj_m
        end
        L = [L,ll];

        % Plot lower bound-----------------------------------------------------
        if verbose        
            set(0,'CurrentFigure',fig_L);                
            plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
            plot(0:numel(L) - 1,L,'b.','markersize',10); hold off                        
        end
    end
        
    % Check convergence----------------------------------------------------
%     [abs(L(end)-L(end - 1)) 2*1e-4*Nm 2*1e-5*Nm]
    if ~(abs(L(end)-L(end - 1))>2*tolmain*Nm)        
        fprintf('==============================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==============================================\n')                        
        break;
    end
end
%==========================================================================

%==========================================================================
function [obj,munum1,muden1] = update_subject_pars(obj,logtpm,iter)

if obj.doaff && ((obj.use_tpm && iter==1) || (~obj.use_tpm && iter==3))
    % Affinely register template to 1st image (n=1)
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,logtpm,obj.Affine,'mni');            
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp, obj.fwhm,      logtpm,obj.Affine,'mni');                                    
    obj.doaff  = 0;
end                 

if ~obj.dodef0
    % Do not estimate deformations
    obj.dodef = 0;
elseif obj.dodef0 && obj.use_tpm && iter==1
    % If a template is used, start estimating deformations from 1st
    % iteration
    obj.dodef = 1;
elseif ~obj.use_tpm && ~obj.dotpm        
    % If no template, and no template update, do not estimate deformations
    obj.dodef = 0;
end

if ~obj.use_tpm && iter==1 && obj.dotpm   
    % 1st iteration
    obj.lkp   = 1:max(obj.lkp); % One Gaussian per tissue
    obj.dodef = 0;              % No deformation update for 1st iteration
elseif ~obj.use_tpm && iter==2 && obj.dotpm                      
    % 2nd iteration
    obj.lkp   = repelem(1:max(obj.lkp),obj.nlkp); % Use nlkp Gaussians per tissue

    if obj.clear_mog_pars
        % Re-estimate cluster parameters based on template constructed after 1st iteration
        if isfield(obj,'mg'), obj = rmfield(obj,'mg'); end
        if isfield(obj,'wp'), obj = rmfield(obj,'wp'); end
        if isfield(obj,'mn'), obj = rmfield(obj,'mn'); end
        if isfield(obj,'vr'), obj = rmfield(obj,'vr'); end
        if isfield(obj,'po'), obj = rmfield(obj,'po'); end
        if isfield(obj,'pr'), obj = rmfield(obj,'pr'); end        
        obj.clear_mog_pars        = 0;
    end
elseif ~obj.use_tpm && iter==4 && obj.dotpm       
    % 4th iteration: start estimating deformations
    obj.dodef = obj.dodef0; 
end         

if isfield(obj,'pth_Twarp')
    % To save memory, load deformation from file
    load(obj.pth_Twarp,'Twarp')
    obj.Twarp = Twarp;
end

% Run segmentation algorithm
[obj,munum1,muden1] = spm_preprocx(obj,logtpm);

% try
%     [obj,munum1,muden1] = spm_preprocx(obj,logtpm);
% catch ME
%     fprintf(ME.message)
%     warning('spm_preprocx')
% 
%     obj.ll = 0;
%     obj.nm = 0;
%     munum1 = single(0);
%     muden1 = single(0);    
% end

if isfield(obj,'pth_Twarp')
    % To save memory, delete deformation from result structure
    Twarp = obj.Twarp;
    save(obj.pth_Twarp,'Twarp')
    obj   = rmfield(obj,'Twarp');
end
%==========================================================================