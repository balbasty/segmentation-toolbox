clear; close all;

%==========================================================================
% TODO
% -MRI+CT
% -Healthy+unhealthy
% -Handle missing data
% -Registration should use affine+diffeo
% -Add MRF
% -Parameterisation: mu-f, matrix exponentials
% -Use a smoothing prior for the template update
% -Test on BRATS data
%==========================================================================

addpath(genpath('./code'))

%==========================================================================
S  = Inf; % Number of subjects
Kb = 15;  % Number of classes (if a template is used, then Kb will be set to the number of classes in that template)

% Define input image data, cell array should contain the following:
% {'pth_to_images','healthy'/'non-healthy',number_of_subjects_to_use,'CT'/'MRI','pth_to_tmpdir'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.
imdir = {'/home/mbrud/Data/Parashkev-CT-aged','healthy',S,'CT','./tmp'};

% Options for input data
load_from_tmpdir = 1; % Run algorithm on images contained in the temporary directory
run_on_2D        = 0; % Run algorithm on 2D data instead of 3D

% Run the algorithm in parallel or not
if 1, runpar = Inf;   % Uses maximum numbers of workers available
else  runpar = 0; end

% Preprocessing options
preproc.realign   = 1; % Realign to MNI space
preproc.crop      = 1; % Remove data outside of head
preproc.denoise   = 1; % Denoise (only done if data is CT)
preproc.create_2D = 1; % Create 2D versions of preprocessed input data

% The distance (mm) between samples (for sub-sampling input data--improves speed)
samp = 1.5;

% Segmentation parameters
nlkp      = 1; % Number of gaussians per tissue
use_vbmog = 1; % Use a variational Bayesian mixture model

% What estimates to perform
dobias = 1; % Bias field
doaff  = 1; % Affine registration
dodef  = 1; % Non-linear registration
dopr   = 1; % Intensity priors
dotpm  = 1; % Template

% Different number of iterations and stopping tolerance of algorithm
nitermain = 100;
niter     = 30;
niter1    = 8;
nsubitmog = 20;
nsubitbf  = 1;
nitdef    = 3;
tol       = 1e-4;

% Regularisation for estimating deformations
rparam = [0 0.001 0.5 0.05 0.2]*0.1;

% Template options
Ptpm     = '';   % Path to existing template (set to '' for estimating a template, or get_spm_TPM for using the default SPM one)
vxtpm    = 1.5;  % Voxel size of template to be estimated
deg      = 2;    % Degree of interpolation when sampling template
tiny     = 1e-4; % Strength of Dirichlet prior used in template construction
fwhm_tpm = 0.25; % Ad hoc smoothing of template

% For debugging
verbose = 1;
figix   = 1;

% Define and create some directories
dirTPM   = './TPM';   % For template and estimated intensity priors
dirTwarp = './Twarp'; % For storing deformation fields
if exist(dirTPM,'dir'),   rmdir(dirTPM,'s');   end; mkdir(dirTPM);
if exist(dirTwarp,'dir'), rmdir(dirTwarp,'s'); end; mkdir(dirTwarp);

%==========================================================================
%% Load and process image data
[V,N,S] = load_and_process_images(imdir,load_from_tmpdir,preproc,runpar,run_on_2D);
% show_images_in_V(V)

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

    figTPM = figure(figix + 5);
    figL   = figure(figix + 6);
    
    distFig; drawnow;
end  
    
obj = cell(S,1);
for s=1:S
    obj{s}.image    = V{s};
    obj{s}.biasfwhm = 60*ones(1,N);
    obj{s}.biasreg  = 1e-3*ones(1,N);
    
    obj{s}.dirTwarp = dirTwarp;
    
    obj{s}.use_vbmog = use_vbmog;
    obj{s}.use_tpm   = use_tpm;
    
    obj{s}.lkp  = repelem(1:Kb,nlkp);
    obj{s}.nlkp = nlkp;
    
    obj{s}.Affine  = eye(4);
    obj{s}.reg     = rparam;
    obj{s}.samp    = samp;
    obj{s}.fwhm    = 0;
    obj{s}.verbose = verbose;
    
    obj{s}.dobias = dobias;
    obj{s}.dodef0 = dodef; 
    obj{s}.dotpm  = dotpm;         
    obj{s}.doaff  = doaff;
    
    obj{s}.niter      = niter;
    obj{s}.niter1     = niter1;
    obj{s}.nsubitmog  = nsubitmog;
    obj{s}.nsubitbf   = nsubitbf;
    obj{s}.nitdef     = nitdef;
    obj{s}.niter_stop = 9;
    
    obj{s}.fig    = fig; 
    
    % Allocate initial deformation field (only for multiple subjects)       
    if S>1
        obj{s}.pthTwarp = alloc_Twarp(obj{s},s);
    end
end

%==========================================================================
%% Run algorithm
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
        pr = update_intensity_prior(obj);
        for s=1:S
            obj{s}.pr = pr;
        end      
        save(fullfile(dirTPM,'pr.mat'),'pr'); % save estimated intensity prior
    end
        
    % Load template
    logtpm = spm_load_logpriors8(Plogtpm,tiny,deg,uniform);
    
    if verbose
        % Visualise template  
        phi = double(identity(logtpm.d));
        b   = spm_sample_logpriors8(logtpm,phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3));
        clear phi

        K1 = floor(sqrt(Kb)); K2 = ceil(Kb/K1); 
        set(0,'CurrentFigure',figTPM);                                        
        for i=1:Kb    
            subplot(K1,K2,i);
            imagesc(b{i}(:,:,floor(logtpm.d(3)/2) + 1)'); axis image xy off; colormap(gray);
        end 
        clear b
        drawnow
    end     
              
    % Iterate over subjects------------------------------------------------
    ll = 0;
    parfor (s=1:S,runpar)    
        fprintf('iter=%d, s=%d\n',iter,s);  
        
        [obj{s},munum1,muden1] = update_model_parameters(obj{s},logtpm,iter);
        
        ll = ll + obj{s}.ll; % Sum up likelihoods over all subjects
        
        if iter==1
           Nm = Nm + obj{s}.nm; % Count total number of voxels in all images
        end
        
        munum = munum + munum1;
        muden = muden + muden1;
    end
           
    L(iter + 1) = ll;

    if verbose
        % Plot lower bound
        set(0,'CurrentFigure',figL);                
        plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
        plot(0:numel(L) - 1,L,'b.','markersize',10); hold off                        
    end
    
    if ~(abs(L(end)-L(end - 1))>2*tol*Nm)
        % Finished
        fprintf('==========================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==========================================\n')                        
        break;
    end
end

%==========================================================================
%% Clean-up estimated template

%==========================================================================