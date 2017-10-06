clear; % close all;

%==========================================================================
% TODO
% -mean correct affine and vecocities (requires different params...)
% -handle missing data
% -registration should use affine+diffeo
% -MRI+CT
% -Healthy+unhealthy
% -MRF?
% -Parameterisation: mu-f, matrix exponentials
% -Test on BRATS data
%==========================================================================

addpath('./util')
addpath('./spm')

% Image data
% imdir  = '/home/mbrud/Data/IXI-subjects/';

% tmpdir = '/home/mbrud/Data/CT-healthy/';
% tmpdir = '/home/mbrud/Desktop/CT-w-lesion/cropped/';
% tmpdir = '/home/mbrud/Desktop/CT-w-lesion/denoised/';
    
% tmpdir = '/home/mbrud/Data/tmpdir-spm_preproc/';

% tmpdir = '/home/mbrud/Desktop/MRI-slices/MRI-slices/';
tmpdir = '/home/mbrud/Desktop/CT-slices/CT-slices-healthy/';
% tmpdir = '/home/mbrud/Desktop/CT-slices/CT-slices-lesion/';

% tmpdir = '/home/mbrud/Desktop/tmp-mri/';
% tmpdir = '/home/mbrud/Desktop/tpm_ixi/IXI012-HH-1211-3DBRAINIXMADisoTFE12_-s3T111_-0301-00003-000001-01.nii';

dirTwarp = './Twarp/';
if exist(dirTwarp,'dir')
    rmdir(dirTwarp,'s');
end
mkdir(dirTwarp);

S  = 16;
Kb = 10;

use_tpm   = 0;
use_vbmog = 1;

nlkp  = 1;
samp  = 2;
vxtpm = 1.5;
deg   = 3; 

dobias = 1;
doaff  = 1;
dodef  = 1;
dopr   = 1;
dotpm  = 1;

rparam     = [0 0.001 0.5 0.05 0.2];
avg_affreg = 1;

tol = 1e-4;

nitermain = 50;
niter     = 30;
niter1    = 8;
nsubitmog = 20;
nsubitbf  = 1;
nitdef    = 3;

tiny = 1e-4;

verbose          = 0;
figix            = 1;
load_from_tmpdir = 1;

runpar = get_nbr_of_cores;
% runpar = 0;

%% Load image data
if load_from_tmpdir
    [V,N,S] = get_V(tmpdir,S);
else
    [V,N,S] = get_V(imdir,S);
    V       = mktmpdir(V,tmpdir);
    V       = reg_and_reslice(V);
end

%% Initialise template
[Plogtpm,Kb] = init_logtpm(use_tpm,V,Kb,vxtpm);

%% Initialise algorithm i/o
lkp = repelem(1:Kb,nlkp);

fig = cell(3,1);
if verbose>1 && ~runpar
%     spm_figure('Create','Interactive');
    
    if verbose>2        
        for i=1:size(fig,1)
           fig{i} = figure((figix - 1)*5 + i);
           clf(fig{i})
        end       
    end
end

if dotpm
    figtpm = figure((figix - 1)*5 + 4);
    clf(figtpm)    
end

figL = figure((figix - 1)*5 + 5);
clf(figL)    

distFig;

obj = cell(S,1);
for s=1:S
    obj{s}.image    = V{s};
    obj{s}.biasfwhm = 60*ones(1,N);
    obj{s}.biasreg  = 1e-3*ones(1,N);
    
    obj{s}.use_vbmog = use_vbmog;
        
    obj{s}.lkp = lkp;
    
    obj{s}.Affine  = eye(4);
    obj{s}.reg     = rparam;
    obj{s}.samp    = samp;
    obj{s}.fwhm    = 0;
    obj{s}.verbose = verbose;

    obj{s}.dobias = dobias;
    if use_tpm
        obj{s}.dodef = dodef;
    else
        obj{s}.dodef = 0;
    end
    obj{s}.dotpm  = dotpm;      
    
    obj{s}.niter     = niter;
    obj{s}.niter1    = niter1;
    obj{s}.nsubitmog = nsubitmog;
    obj{s}.nsubitbf  = nsubitbf;
    obj{s}.nitdef    = nitdef;
    
    if verbose>2
        obj{s}.fig = fig;        
    end
end

%% Run algorithm
Nm = 0;
L  = -Inf;
for iter=1:nitermain
    
    % Update the template--------------------------------------------------
    if dotpm && iter>1        
        Nii = nifti(Plogtpm);
        dm  = size(Nii.dat(:,:,:));
        if numel(dm)==2
            dm(3) = 1;
        end
        
        fwhm          = samp/vxtpm/2;
        [munum,muden] = smooth_template(munum,muden,dm,fwhm);

        nlogtpm = log(munum./muden + tiny);
        clear munum muden
        
        mu = safe_softmax(nlogtpm');
        mu = reshape(mu,[dm Kb]);        
        z  = floor(dm(3)/2 + 1);
        K1 = floor(sqrt(Kb));
        K2 = ceil(Kb/K1); 
        set(0,'CurrentFigure',figtpm);              
        for i=1:Kb
            subplot(K1,K2,i);
            imagesc(mu(:,:,z,i)'); axis image xy off; title(['k=' num2str(i)]); colormap(gray);
        end  
        drawnow
        clear mu
                
        Nii.dat(:,:,:,:) = reshape(nlogtpm',[dm Kb]);
        clear Nii nlogtpm
    end
    
    % Update intensity "priors"--------------------------------------------
    if dopr && use_vbmog && iter>2
        pr = update_pr(obj);
        for s=1:S
            obj{s}.pr = pr;
        end      
    end
        
    % Load template
    logtpm = spm_load_logpriors8(Plogtpm,tiny,deg);
    
    % Affine registration of template to subject
    if doaff && (use_tpm || (~use_tpm && iter==2))
        parfor (s=1:S,runpar)   
%         for s=1:S
            obj{s}.Affine = eye(4);
            obj{s}.Affine = spm_logmaff8(obj{s}.image(1),obj{s}.samp,(obj{s}.fwhm+1)*16,logtpm,obj{s}.Affine,'mni');            
            obj{s}.Affine = spm_logmaff8(obj{s}.image(1),obj{s}.samp, obj{s}.fwhm,      logtpm,obj{s}.Affine,'mni');                                    
        end
        
        % Mean correct affine transform parameters
        if S>1 && avg_affreg
            r = zeros([12 1 S]);
            for s=1:S
                r(:,:,s) = M2P(obj{s}.Affine);
            end
            ravg = mean(r,3);        

            for s=1:S
                obj{s}.Affine = P2M(r(:,:,s) - ravg);
            end
            clear r
        end  
    end        
    
    % For updating template
    if dotpm
        munum = single(0); muden = single(0);            
    end
              
    % Iterate over subjects------------------------------------------------
    ll = 0;
    parfor (s=1:S,runpar)    
%     for s=1:S    
        fprintf('iter=%d, s=%d\n',iter,s);  
        
        % Allocate initial deformation field (only for multiple subjects)       
        if (use_tpm && S>1 && iter==1) || (dotpm && S>1 && iter>=2)
            d0 = obj{s}.image(1).dim;
            vx = sqrt(sum(obj{s}.image(1).mat(1:3,1:3).^2));
            sk = max([1 1 1],round(obj{s}.samp*[1 1 1]./vx));
            x0 = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
            z0 = 1:sk(3):d0(3);
            d  = [size(x0) length(z0)];
            d0 = []; vx = []; sk = []; x0 = []; z0 = [];
            
            Twarp           = zeros([d 3],'single');    
            pthTwarp        = fullfile(dirTwarp,['Twarp' num2str(s) '.nii']);     
            create_nii(pthTwarp,Twarp,eye(4),'float32','Velocity field')
            obj{s}.pthTwarp = pthTwarp;
            Twarp           = [];
        end                

        if ~use_tpm && iter==1 
            obj{s}.lkp = 1:Kb;
            
            [nm,wp,mn,vr] = spm_InitGaussians(obj{s}.image,Kb,obj{s}.samp);            
            Nm            = Nm + nm;
            
            if obj{s}.use_vbmog
                m0 = zeros(N,Kb);
                b0 = ones(1,Kb);
                W0 = zeros(N,N,Kb);
                for k=1:Kb
                    W0(:,:,k) = 200*eye(N,N);
                end
                n0 = 20*ones(1,Kb);

                % Use 'responsibilities' from initialization to set sufficient statistics
                mm0 = nm*wp;
                mm1 = mn;
                mm2 = vr;

                b = b0 + mm0;
                n = n0 + mm0;
                m = zeros(N,Kb);
                W = zeros(N,N,Kb);
                for k=1:Kb
                    m(:,k)   = (b0(k)*m0(:,k) + mm0(k).*mm1(:,k))./b(k);

                    invW0    = inv(W0(:,:,k));
                    mlt1     = b0(k).*mm0(k)/(b0(k) + mm0(k));
                    diff1    = mm1(:,k) - m0(:,k);
                    W(:,:,k) = inv(invW0 + mm0(k)*mm2(:,:,k) + mlt1*(diff1*diff1'));
                end  

                obj{s}.pr.m = m0;
                obj{s}.pr.b = b0;
                obj{s}.pr.n = n0;
                obj{s}.pr.W = W0;
                
                obj{s}.po.m = m;
                obj{s}.po.b = b;
                obj{s}.po.n = n;
                obj{s}.po.W = W;
                
                obj{s}.mg = ones(Kb,1);
            else               
                obj{s}.mn = mn;
                obj{s}.vr = vr;    
                
                obj{s}.mg = ones(Kb,1);
            end
        elseif ~use_tpm && iter==2 && dotpm                      
            if isfield(obj{s},'mg'), obj{s} = rmfield(obj{s},'mg'); end
            if isfield(obj{s},'wp'), obj{s} = rmfield(obj{s},'wp'); end
            if isfield(obj{s},'mn'), obj{s} = rmfield(obj{s},'mn'); end
            if isfield(obj{s},'vr'), obj{s} = rmfield(obj{s},'vr'); end
            if isfield(obj{s},'po'), obj{s} = rmfield(obj{s},'po'); end
            if isfield(obj{s},'pr'), obj{s} = rmfield(obj{s},'pr'); end        
            if isfield(obj{s},'Tbias'), obj{s} = rmfield(obj{s},'Tbias'); end 
            if isfield(obj{s},'Twarp'), obj{s} = rmfield(obj{s},'Twarp'); end 
        
            obj{s}.lkp   = repelem(1:Kb,nlkp);           
            obj{s}.dodef = dodef;
        elseif iter==1
            nm = spm_InitGaussians(obj{s}.image,Kb,obj{s}.samp);            
            Nm = Nm + nm;
        end         
        
        [obj{s},munum1,muden1] = spm_preprocx(obj{s},logtpm);
        ll                     = ll + obj{s}.ll;
        
        if dotpm
            munum = munum + munum1;
            muden = muden + muden1;
        end
    end
           
    L(iter + 1) = ll;

    set(0,'CurrentFigure',figL);                
    plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
    plot(0:numel(L) - 1,L,'b.','markersize',10); hold off                        
    
    if ~(abs(L(end)-L(end - 1))>2*tol*Nm) && iter>20
        % Finished
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        break;
    end
end
