clear; close all;

addpath(genpath('./code'))

%==========================================================================
S  = 64; % Number of subjects
Kb = 15; % Number of classes (if a template is used, then Kb will be set to the number of classes in that template)

%--------------------------------------------------------------------------
% Define input image data, cell array should contain the following:
% {'pth_to_images','healthy'/'non-healthy',number_of_subjects_to_use,'CT'/'MRI','pth_to_tmpdir'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.
% imobj = {'/home/mbrud/Data/Slices/CT-hemorrhage-2D/',S,'CT','healthy'};

imobj{1} = {'/home/mbrud/Data/Slices/IXI-2D/',          S,'MRI','healthy'};
imobj{2} = {'/home/mbrud/Data/Slices/CT-hemorrhage-2D/',S,'CT', 'healthy'};

%--------------------------------------------------------------------------
% Run the algorithm in parallel or not
if 1, runpar = Inf;   % Uses maximum numbers of workers available
else  runpar = 0; end

%--------------------------------------------------------------------------
% Preprocessing options
preproc.do_preproc = 0; % Do preprocessing on input images
preproc.realign    = 1; % Realign to MNI space
preproc.crop       = 1; % Remove data outside of head
preproc.denoise    = 1; % Denoise (only done if data is CT)
preproc.is_DICOM   = 0; % If input images are DICOM, converts DICOM to Nifti % TODO (work in progress)

%--------------------------------------------------------------------------
% The distance (mm) between samples (for sub-sampling input data--improves speed)
samp = 2;

%--------------------------------------------------------------------------
% Segmentation parameters
nlkp      = 1; % Number of gaussians per tissue
use_vbmog = 1; % Use a variational Bayesian mixture model

%--------------------------------------------------------------------------
% What estimates to perform
dobias = 1; % Bias field
doaff  = 1; % Affine registration
dodef  = 1; % Non-linear registration
dopr   = 1; % Intensity priors
dotpm  = 1; % Template

%--------------------------------------------------------------------------
% Different number of iterations and stopping tolerance of algorithm
nitermain = 100;
niter     = 30;
niter1    = 8;
nsubitmog = 20;
nsubitbf  = 1;
nitdef    = 3;
tol       = 1e-4;

%--------------------------------------------------------------------------
% Regularisation for estimating deformations
rparam = [0 0.001 0.5 0.05 0.2]*0.1;

%--------------------------------------------------------------------------
% Template options
Ptpm     = '';   % Path to existing template (set to '' for estimating a template, or get_spm_TPM for using the default SPM one)
vxtpm    = 1.5;  % Voxel size of template to be estimated
deg      = 2;    % Degree of interpolation when sampling template
tiny     = 1e-4; % Strength of Dirichlet prior used in template construction
fwhm_tpm = 0.25; % Ad hoc smoothing of template

%--------------------------------------------------------------------------
% For debugging
verbose = 1;
figix   = 1;

%--------------------------------------------------------------------------
% Define and create some directories
dirTPM   = './TPM';   % For template and estimated intensity priors
dirTwarp = './Twarp'; % For storing deformation fields
if exist(dirTPM,'dir'),   rmdir(dirTPM,'s');   end; mkdir(dirTPM);
if exist(dirTwarp,'dir'), rmdir(dirTwarp,'s'); end; mkdir(dirTwarp);

%==========================================================================
%% Load and process image data
[V,M,S,N] = load_and_process_images(imobj,preproc,runpar);

%==========================================================================
%% Initialise template
[Plogtpm,Kb,uniform,use_tpm] = init_template(Ptpm,V,Kb,vxtpm,dirTPM);

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
    
obj = cell(1,M);
for m=1:M    

    obj{m} = cell(1,S(m));
    for s=1:S(m)
        obj{m}{s}.image    = V{m}{s};
        obj{m}{s}.biasfwhm = 60*ones(1,N(m));
        obj{m}{s}.biasreg  = 1e-3*ones(1,N(m));

        obj{m}{s}.dirTwarp = dirTwarp;

        obj{m}{s}.use_vbmog = use_vbmog;
        obj{m}{s}.use_tpm   = use_tpm;

        obj{m}{s}.lkp  = repelem(1:Kb,nlkp);
        obj{m}{s}.nlkp = nlkp;

        obj{m}{s}.Affine  = eye(4);
        obj{m}{s}.reg     = rparam;
        obj{m}{s}.samp    = samp;
        obj{m}{s}.fwhm    = 0;
        obj{m}{s}.verbose = verbose;

        obj{m}{s}.dobias = dobias;
        obj{m}{s}.dodef0 = dodef; 
        obj{m}{s}.dotpm  = dotpm;         
        obj{m}{s}.doaff  = doaff;

        obj{m}{s}.niter      = niter;
        obj{m}{s}.niter1     = niter1;
        obj{m}{s}.nsubitmog  = nsubitmog;
        obj{m}{s}.nsubitbf   = nsubitbf;
        obj{m}{s}.nitdef     = nitdef;
        obj{m}{s}.niter_stop = 9;

        obj{m}{s}.descrip = imobj{m}{3};
        obj{m}{s}.healthy = imobj{m}{4}; 

        obj{m}{s}.fig = fig; 
              
        if S(m)>1
            % Allocate initial deformation field (only for multiple subjects) 
            obj{m}{s}.pthTwarp = alloc_Twarp(obj{m}{s},m,s);
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
L  = -Inf; % Lower bound of model
for iter=1:nitermain
    
    % Update the template--------------------------------------------------                     
    if dotpm && iter>1   
        uniform = false;
        
        Nii = nifti(Plogtpm);
        dm  = size(Nii.dat(:,:,:));
        if numel(dm)==2
            dm(3) = 1;
        end
        
        % Smooth template
        [munum,muden] = smooth_template(munum,muden,dm,fwhm_tpm);

        nlogtpm = log(munum./muden + tiny);
        clear munum muden
               
        Nii.dat(:,:,:,:) = reshape(nlogtpm',[dm Kb]);
        clear Nii nlogtpm        
    end
    
    munum = single(0); muden = single(0);   
    
    % Update intensity prior-----------------------------------------------
    if dopr && use_vbmog && iter>=3
        for m=1:M
            pr = update_intensity_prior(obj{m});
            for s=1:S(m)
                obj{m}{s}.pr = pr;
            end      
            save(fullfile(dirTPM,['pr_m' num2str(m) '.mat']),'pr'); 
        end
    end
        
    % Load template
    logtpm = spm_load_logpriors8(Plogtpm,tiny,deg,uniform);
    
    if verbose
        % Visualise template  
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
              
    % Iterate over subjects------------------------------------------------
    ll = 0;
    for m=1:M
        obj_m = obj{m};

%         for s=1:S
        parfor (s=1:S(m),runpar)    
            fprintf('iter=%d, m=%d, s=%d\n',iter,m,s);  

            [obj_m{s},munum1,muden1] = update_model_parameters(obj_m{s},logtpm,iter);

            ll = ll + obj_m{s}.ll; % Sum up likelihoods over all subjects

            if iter==1
               Nm = Nm + obj_m{s}.nm; % Count total number of voxels in all images
            end

            munum = munum + munum1;
            muden = muden + muden1;
        end

        obj{m} = obj_m;
        clear obj_m
    end
    L(iter + 1) = ll;

    if verbose
        % Plot lower bound
        set(0,'CurrentFigure',fig_L);                
        plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
        plot(0:numel(L) - 1,L,'b.','markersize',10); hold off                        
    end
    
    if ~(abs(L(end)-L(end - 1))>2*tol*Nm)
        % Finished
        fprintf('==============================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==============================================\n')                        
        break;
    end
end

%==========================================================================
%% Clean-up estimated template

%==========================================================================