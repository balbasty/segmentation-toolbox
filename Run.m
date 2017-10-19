clear; close all;

%==========================================================================
% TODO
% -MRI+CT (need CT and MRI of the same subject)
% -Healthy+unhealthy
% -handle missing data
% -registration should use affine+diffeo
% -Add MRF
% -Parameterisation: mu-f, matrix exponentials
% -Test on BRATS data
%==========================================================================

addpath('./util')
addpath('./spm')
addpath('./preproc')

% Image data
load_from_tmpdir = 0;

S = 1000; % number of subjects

imdir  = {};
tmpdir = {};

imdir{end + 1}  = {'/home/mbrud/Desktop/Parashkev-CT-aged/','healthy',S,'CT'};
% imdir{end + 1}  = {'/home/mbrud/Desktop/CT-slices/CT-slices-lesion/','healthy',S,'CT'};
% imdir{end + 1}  = {'/home/mbrud/Desktop/nifti-toolbox/resize-imgs/resized-data/','healthy',S,'CT'};
tmpdir{end + 1} = {'./tmp/','healthy',S,'CT'};
% tmpdir{end + 1} = {'./tmp_2D/','healthy',S,'CT'};
% tmpdir{end + 1} = {imdir{end}{1},'healthy',S,'CT'};

% Algorithm parameters
Kb    = 15;
Ksemi = 0;
nlkp  = 1;

use_tpm   = 0;
use_vbmog = 1;

samp  = 1.5;
vxtpm = 1.5;
deg   = 3; 

preproc.realign   = 1;
preproc.crop      = 1;
preproc.denoise   = 1;
preproc.create_2D = 1;

dobias = 1;
doaff  = 1;
dodef  = 1;
dopr   = 1;
dotpm  = 1;

rparam = [0 0.001 0.5 0.05 0.2]*0.1;

fwhm_tpm = 0.5; % ad hoc...

tol = 1e-4;

nitermain = 50;
niter     = 30;
niter1    = 8;
nsubitmog = 20;
nsubitbf  = 1;
nitdef    = 3;

tiny = 1e-4;

verbose = 1;
figix   = 1;

if 1
    runpar = get_nbr_of_cores;
else
    runpar = 0;
end

% Make directories
dirTPM    = './TPM/';
dirTwarp  = './Twarp/';

if exist(dirTPM,'dir')
    rmdir(dirTPM,'s');
end
mkdir(dirTPM);

if exist(dirTwarp,'dir')
    rmdir(dirTwarp,'s');
end
mkdir(dirTwarp);

%% Load and process image data

% Assumes that images are organised into per subject subfolders, i.e. if
% subject s=1 has a T1 and a T2 image, then they are in a subfolder S1.
[V,N,S] = load_and_process_images(imdir,tmpdir,load_from_tmpdir,preproc);
% show_images_in_V(V)

% Determine if a semi-supervised approach should be used, based on knowing if
% input images contain pathology or not
healthy_nonhealthy = {};
for i=1:numel(imdir)
    healthy_nonhealthy{i} = imdir{i}{2};
end
s1 = strcmp('healthy',healthy_nonhealthy);
s2 = strcmp('nonhealthy',healthy_nonhealthy);

semi_healthy = ~(all(s1 == s1(1)) && all(s2 == s2(1)));

for s=1:S
    for n=1:N
        V{s}(n).brain_is = imdir{i}{2};
        V{s}(n).descrip  = imdir{i}{4};
    end    
end

%% Initialise template
[Plogtpm,Kb,uniform] = init_logtpm(use_tpm,V,Kb,vxtpm,dirTPM);
lkp                  = repelem(1:Kb,nlkp);

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
    
    obj{s}.lkp  = lkp;
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
    
    % Determine if a semi-supervised approach should be used, based on knowing if
    % input images contain pathology or not
    if semi_healthy && strcmp(obj{s}.image.brain_is,'healthy') && Ksemi
        lkpsemi        = fliplr(1:Kb);
        obj{s}.lkpsemi = lkpsemi(1:Ksemi);
    else
        obj{s}.lkpsemi = [];
    end

    obj{s}.runpar = runpar; 
    obj{s}.fig    = fig; 
    
    % Allocate initial deformation field (only for multiple subjects)       
    if S>1
        obj{s}.pthTwarp = alloc_Twarp(obj{s},s);
    end
end

%% Run algorithm

Nm = 0; % Total # voxels in all images
L  = -Inf;
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
        
        ll = ll + obj{s}.ll;
        
        if iter==1
           Nm = Nm + obj{s}.nm; % Count total number of voxels in all images
        end
        
        munum = munum + munum1;
        muden = muden + muden1;
    end
           
    L(iter + 1) = ll;

    if verbose
        set(0,'CurrentFigure',figL);                
        plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
        plot(0:numel(L) - 1,L,'b.','markersize',10); hold off                        
    end
    
    if ~(abs(L(end)-L(end - 1))>2*tol*Nm)
        % Finished
        fprintf('==========================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==========================================\n')                        
%         break;
    end
end

%% Clean-up estimated template