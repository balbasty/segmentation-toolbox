function spm_preprocx_run(obj,im,K)

%==========================================================================
% Make some directories
%==========================================================================

if obj.run_on_holly
    [~,obj.dir_data,~] = read_directory_details('directory_details.txt',obj.holly_jnam);    
else
    obj.holly_jnam     = '';
end

if ~exist(obj.dir_data,'dir'), mkdir(obj.dir_data); end
obj.dir_res = fullfile(obj.dir_data,obj.dir_res); 
if exist(obj.dir_res,'dir'), rmdir(obj.dir_res,'s'); end; mkdir(obj.dir_res);

%==========================================================================
% Load and process image data
%==========================================================================

[V,labels] = load_and_process_images(obj,im);

%==========================================================================
% Initialise template
%==========================================================================

[obj,K] = init_tpm(obj,V,K);

%==========================================================================
% Initialise some debugging output
%==========================================================================

[fig_seg,fig_tpm,fig_L,obj.verbose] = create_fig(obj);

%==========================================================================
% Initialise algorithm i/o 
%==========================================================================

nitermain    = obj.nitermain;
tolmain      = obj.tolmain;
dotpm        = obj.dotpm;
vb           = obj.vb;
use_mog      = ~isempty(obj.lkp);
dopr         = obj.dopr && use_mog && vb;
pth_logtpm   = obj.pth_logtpm;
pr_dirichlet = obj.pr_dirichlet;
deg          = obj.deg;
fwhmtpm      = obj.fwhmtpm;
dir_data     = obj.dir_data;
dir_res      = obj.dir_res;
run_on_holly = obj.run_on_holly;
num_workers  = obj.num_workers;
holly_jnam   = obj.holly_jnam;
crop_bb      = obj.crop_bb;

pth_obj = get_pth_obj(obj,V,K,im,labels,run_on_holly);
clear V obj

%==========================================================================
% Prepare for running the algorithm on the Holly cluster
% OBS! For this to work, this code needs to be located on the Holly server
%==========================================================================
   
if run_on_holly
    holly = init_holly(pth_obj,dir_data,holly_jnam);
else
    holly = [];
end

%==========================================================================
% Run algorithm
%==========================================================================

print_algorithm_progress('started');

L = -Inf; % Lower bound of complete model
for iter=1:nitermain
    fprintf('iter=%d======================\n',iter);  
                   
    % Update subject specific parameters (clusters, bias, affine, deformations, template derivatives)
    [L,tpmnum,tpmden,Nm,holly] = update_subjects(pth_obj,L,num_workers,run_on_holly,holly,fig_seg,pth_logtpm); 
     
    % Plot objective value
    plot_objval(L,fig_L); 
        
    if dotpm
        % Update template 
        update_global_tpm(tpmnum,tpmden,pth_logtpm,pr_dirichlet,fwhmtpm,pth_obj,crop_bb,iter);                
              
        show_tpm(fig_tpm,pth_logtpm,deg);
    end  
    
    if dopr
        % Update intensity prior
        update_global_prior(pth_obj,im);
    end
        
    % Check convergence
    fprintf('L=%0.0f\n',L(end));    
    if ~(abs(L(end)-L(end - 1))>2*tolmain*Nm) && iter>9
        fprintf('==============================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==============================================\n')                        
          
        % Convert log template to probability template and save to disk
        logtpm2tpm(pth_logtpm,dir_res);
    
        break;
    end
end

print_algorithm_progress('finished');
%==========================================================================

%==========================================================================
function [L,tpmnum,tpmden,Nm,holly] = update_subjects(pth_obj,L,num_workers,run_on_holly,holly,fig,pth_logtpm)
V      = spm_vol(pth_logtpm);
dm_tpm = V(1).dim;
K      = numel(V);
clear V

if run_on_holly
    % Run subject specific jobs on the FIL cluster (Holly)
    %----------------------------------------------------------------------            
    holly = parfor_holly(holly);
else 
    % Run subject specific jobs using MATLAB parfor
    %----------------------------------------------------------------------
    parfor_matlab(pth_obj,num_workers,fig);
end

% Read results from estimations
%--------------------------------------------------------------------------
tic;
M      = numel(pth_obj);
tpmnum = zeros([dm_tpm K],'single'); tpmden = tpmnum; ll = 0; Nm = 0; totS = zeros(1,M); sDC = cell(1,M); sDC(1,:) = {0};
for m=1:M
    S = numel(pth_obj{m});             
    for s=1:S
        obj = load(pth_obj{m}{s},'-mat'); 

        if obj.status==0   
            if obj.dotpm
                [~,ix_x,ix_y,ix_z] = bb_info(obj.bb,dm_tpm);

                if run_on_holly
                    Nii = nifti(obj.pth_tpmnum_l);
                else
                    Nii = nifti(obj.pth_tpmnum);
                end                
                tpmnum1 = Nii.dat(:,:,:,:); 
                tpmnum(ix_x,ix_y,ix_z,:) = tpmnum(ix_x,ix_y,ix_z,:) + tpmnum1;

                if run_on_holly
                    Nii = nifti(obj.pth_tpmden_l);
                else
                    Nii = nifti(obj.pth_tpmden);
                end
                tpmden1 = Nii.dat(:,:,:,:); 
                tpmden(ix_x,ix_y,ix_z,:) = tpmden(ix_x,ix_y,ix_z,:) + tpmden1; 
            end
            
            sDC{m}  = sDC{m} + obj.DC;
            totS(m) = totS(m) + 1;
        else
            fprintf(2,['Error for image: ' obj.image(1).fname '\n']);
        end

        ll = ll + obj.ll;
        Nm = Nm + obj.nm;
    end
end
L = [L,ll];

% For setting the DC component of all the bias fields so that they
% average to 0 (used for global signal normalisation)
%--------------------------------------------------------------------------
avgDC = cell(1,M);
for m=1:M   
    avgDC{m} = sDC{m}/totS(m);
    
    if 1
       fprintf('avgDC{%d} = %s\n',m,sprintf('%0.2f ',avgDC{m}));  
    end
end

for m=1:M
    S = numel(pth_obj{m});             
    for s=1:S
        obj    = load(pth_obj{m}{s},'-mat'); 
        obj.DC = avgDC{m};
        save(pth_obj{m}{s},'-struct','obj')                    
    end
end

fprintf('Elapsed time (update_subjects): %0.1f s\n',toc);  
%==========================================================================

%==========================================================================
function parfor_matlab(pth_obj,num_workers,fig)
% Estimate on all subjects
M = numel(pth_obj);
tic;
manage_parpool(num_workers);
for m=1:M      
    fprintf('Updating parameters for subject group m=%d\n',m);  
    
    pth_obj_m = pth_obj{m};
    S         = numel(pth_obj_m);    
    parfor (s=1:S,num_workers)                                  
%     for s=1:S
        obj = load(pth_obj_m{s});            
        
        % Run segmentation routine         
        obj = segment(obj,fig);                                       

        save_in_parfor(pth_obj_m{s},obj,'-struct');
    end     
end    
fprintf('Elapsed time (parfor_matlab): %0.1f s\n',toc);  
%==========================================================================

%==========================================================================
function update_global_prior(pth_obj,im,verbose)
if nargin<3, verbose = false; end

tic
M = numel(im);

% Map modalities
% img_mod = cell(1,M);
% for m=1:M
%     img_mod{m} = im{m}{3};
% end
% [~,M1,img_mod] = unique(img_mod);
% M1             = numel(M1);
img_mod = 1:M;
M1      = M;

% Load posteriors and priors from all subjects and store in struct
obj = cell(1,M1);
for m=1:M
    m1     = img_mod(m);
    S      = numel(pth_obj{m});    
    cnt    = numel(obj{m1});
    for s=1:S
        obj1 = load(pth_obj{m}{s},'-mat','mog','status');
        if obj1.status==0
            cnt             = cnt + 1;
            obj{m1}{cnt}.po = obj1.mog.po;
            obj{m1}{cnt}.pr = obj1.mog.pr;            
        end
    end
end

% Update prior based on posteriors and previous priors, then save new prior
fprintf('Updating Gaussian-Wishart hyperparameters...\n')
pr = cell(1,M1);
for m=1:M1
    if ~isempty(obj{m})
        pr{m} = update_intensity_prior(obj{m});   
          
        if verbose
            N = size(pr{m}.m,1);
            for n=1:N
               fprintf('pr(m=%d).m%d = [%0.2f, %s%0.2f]\n',m,n,pr{m}.m(n,1),sprintf('%0.2f, ',pr{m}.m(n,2:end-1)),pr{m}.m(n,end)); 
            end
            fprintf('pr(m=%d).b = [%0.2f, %s%0.2f]\n',m,pr{m}.b(1),sprintf('%0.2f, ',pr{m}.b(2:end-1)),pr{m}.b(end)); 
            fprintf('pr(m=%d).n = [%0.2f, %s%0.2f]\n',m,pr{m}.n(1),sprintf('%0.2f, ',pr{m}.n(2:end-1)),pr{m}.n(end)); 
        end
    end
end
clear obj

for m=1:M  
    m1        = img_mod(m);
    S         = numel(pth_obj{m});
    for s=1:S     
        obj        = load(pth_obj{m}{s});  
        obj.mog.pr = pr{m1};   
        save(pth_obj{m}{s},'-struct','obj')
    end              
end
fprintf('Elapsed time (update_global_prior): %0.1f s\n',toc);
%==========================================================================        

%==========================================================================
function update_global_tpm(tpmnum,tpmden,pth_logtpm,pr_dirichlet,fwhmtpm,pth_obj,crop_bb,iter,softmax_tpm,dir_res)
if nargin<9, softmax_tpm = false; end

fprintf('Updating TPMs...\n')

% Smooth template
[tpmnum,tpmden] = smooth_tpm(tpmnum,tpmden,fwhmtpm);

% Update template
logtpm = log(tpmnum./tpmden + pr_dirichlet);
clear tpmnum tpmden

% Save updated template
Nii = nifti(pth_logtpm);
Nii.dat(:,:,:,:) = logtpm;        
clear Nii

if crop_bb && iter>=3 && size(logtpm,3)>1
    % Crop template according to subject bounding boxes
    M = numel(pth_obj);
    S = 0;
    for m=1:M, S = S + numel(pth_obj{m}); end

    bb  = [];
    cnt = 1;
    for m=1:M
        S = numel(pth_obj{m});    
        for s=1:S
            obj = load(pth_obj{m}{s},'-mat','bb','status');
            if obj.status==0 && ~isempty(obj.bb)
                bb1         = obj.bb;
                bb1         = [bb1(3) bb1(4);bb1(1) bb1(2);bb1(5) bb1(6)];
                bb(:,:,cnt) = bb1;
                cnt         = cnt + 1;
            end
        end
    end

    mn_bb    = min(bb,[],3);
    mx_bb    = max(bb,[],3);
%     nbb      = [mn_bb(:,1) mx_bb(:,2)];
    nbb      = [mx_bb(:,1) mn_bb(:,2)];
    nbb(3,1) = min(0,nbb(3,1)); % To not accidentally remove the neck

    V0  = spm_vol(pth_logtpm);
    od = V0(1).dim;

    for k=1:numel(V0)
        subvol(V0(k),nbb','tmp');        
    end

    delete(pth_logtpm);
    [pth,nam,ext] = fileparts(V0(1).fname);
    fname         = fullfile(pth,['tmp' nam ext]);
    movefile(fname,pth_logtpm);

    V  = spm_vol(pth_logtpm);
    nd = V(1).dim;

    fprintf('size(tpm)=[%d %d %d] => size(ntpm)=[%d %d %d]\n',od(1),od(2),od(3),nd(1),nd(2),nd(3));
end

if softmax_tpm
    % Convert log template to probability template and save to disk
    logtpm2tpm(pth_logtpm,dir_res);
end
%==========================================================================

%==========================================================================
function show_tpm(fig_tpm,pth_logtpm,deg)
if ~isempty(fig_tpm) 
    logtpm = spm_load_logpriors8(pth_logtpm,deg);

    phi = double(identity(logtpm.d));
    b   = spm_sample_logpriors8(logtpm,phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3));
    clear phi

    K  = numel(b);
    K1 = floor(sqrt(K)); 
    K2 = ceil(K/K1); 
    set(0,'CurrentFigure',fig_tpm);                                        
    for i=1:K    
        subplot(K1,K2,i);
        imagesc(b{i}(:,:,floor(logtpm.d(3)/2) + 1)'); axis image xy off; colormap(gray);
        title(['tpm, k=' num2str(i)]);
    end 
    clear b
    drawnow
end
%==========================================================================

%==========================================================================
function plot_objval(L,fig_L)
if ~isempty(fig_L) && numel(L)>=3
    set(0,'CurrentFigure',fig_L);                
    plot(0:numel(L(3:end)) - 1,L(3:end),'b-','LineWidth',1);   hold on            
    plot(0:numel(L(3:end)) - 1,L(3:end),'b.','markersize',10); hold off  
    title('Lower bound')
end
%==========================================================================

%==========================================================================
function pth_obj = get_pth_obj(obj0,V,K,im,labels,run_on_holly)
dir_data = obj0.dir_data;

dir_obj  = fullfile(dir_data,'obj');
if exist(dir_obj,'dir'), rmdir(dir_obj,'s'); end; mkdir(dir_obj);

dir_def  = fullfile(dir_data,'def');
if exist(dir_def,'dir'), rmdir(dir_def,'s'); end; mkdir(dir_def);

dir_tpmnum  = fullfile(dir_data,'tpmnum');
if exist(dir_tpmnum,'dir'), rmdir(dir_tpmnum,'s'); end; mkdir(dir_tpmnum);

dir_tpmden  = fullfile(dir_data,'tpmden');
if exist(dir_tpmden,'dir'), rmdir(dir_tpmden,'s'); end; mkdir(dir_tpmden);

M       = numel(V);
pth_obj = cell(1,M);
for m=1:M       
    S          = numel(V{m});
    pth_obj{m} = cell(1,S);
    for s=1:S  
        obj        = struct;
        obj.status = false;
        
        N         = numel(V{m}{s});
        obj.image = V{m}{s};
        if run_on_holly
            for n=1:N
                fname              = obj.image(n).fname;
                nfname             = ['/' fname(7:end)];
                obj.image(n).fname = nfname;            
            end
        end        
        
        obj.biasfwhm = obj0.biasfwhm*ones(1,N);
        obj.biasreg  = obj0.biasreg*ones(1,N);       
        obj.DC       = zeros(1,N);
        
        obj.use_tpm = obj0.use_tpm;
        
        obj.missing_data = obj0.missing_data;
        
        obj.mog.vb     = obj0.vb;        
        obj.wp_reg     = obj0.wp_reg;
        if obj0.dotpm
            % One Gaussian per tissue for template construction
            obj.lkp = 1:K; 
        else
            if isempty(obj0.lkp)
                obj.lkp = 0;
            elseif isscalar(obj0.lkp)
                obj.lkp = repelem(1:K,obj0.lkp);
            else
                obj.lkp = obj0.lkp;
            end
        end        

        obj.Affine  = eye(4);
%         obj.r       = zeros([12,1]); % Matrix exponentials parameterisation of an affine transformation
               
        obj.reg     = obj0.rparam;
        obj.samp    = obj0.samp;
        obj.fwhm    = 0;
        obj.verbose = obj0.verbose;

        obj.dobias   = obj0.dobias;
        obj.dodef    = obj0.dodef; 
        obj.dodef0   = obj0.dodef; 
        obj.dotpm    = obj0.dotpm;  
        obj.doaff    = obj0.doaff;
        obj.dowp     = obj0.dowp;
        obj.dowp0    = obj0.dowp;
                   
        obj.tolseg    = obj0.tolseg;
        obj.niter     = obj0.niter;
        obj.niter1    = obj0.niter1;
        obj.nsubitmog = obj0.nsubitmog;
        obj.nsubitbf  = obj0.nsubitbf;
        obj.nitdef    = obj0.nitdef;

        obj.descrip = im{m}{3};
        obj.healthy = im{m}{4};         
        
        obj.def_done = false;
        
        obj.ll = 0;
        obj.nm = 0;
              
        obj.bb = [];
        
        obj.pth_logtpm = obj0.pth_logtpm;
        if run_on_holly
            fname          = obj.pth_logtpm;
            nfname         = ['/' fname(7:end)];
            obj.pth_logtpm = nfname;            
        end
           
        obj.deg  = obj0.deg;        
        obj.iter = 0;
               
        % Path to deformations
        pth_def     = fullfile(dir_def,['def-m' num2str(m) '-s' num2str(s) '.nii']);                       
        obj.pth_def = pth_def;        
        if run_on_holly        
            fname       = pth_def;
            nfname      = ['/' fname(7:end)];
            obj.pth_def = nfname;                        
        end    
        
        % Path to template updates
        pth_tpmnum     = fullfile(dir_tpmnum,['tpmnum-m' num2str(m) '-s' num2str(s) '.nii']);  
        obj.pth_tpmnum = pth_tpmnum;                
        if run_on_holly        
            obj.pth_tpmnum_l = pth_tpmnum;
            fname           = pth_tpmnum;
            nfname          = ['/' fname(7:end)];
            obj.pth_tpmnum   = nfname;                        
        end   
        
        pth_tpmden     = fullfile(dir_tpmden,['tpmden-m' num2str(m) '-s' num2str(s) '.nii']); 
        obj.pth_tpmden = pth_tpmden;                
        if run_on_holly        
            obj.pth_tpmden_l = pth_tpmden;
            fname           = pth_tpmden;
            nfname          = ['/' fname(7:end)];
            obj.pth_tpmden   = nfname;                        
        end   
        
        % Allocate labels
%         obj.labels  = labels{m}{s}; 
        
        % Store path to obj        
        pth_obj{m}{s} = fullfile(dir_obj,['obj-m' num2str(m) '-s' num2str(s) '.mat']);
        
        save(pth_obj{m}{s},'-struct','obj')
    end
end
%==========================================================================

%==========================================================================
function [fig_seg,fig_tpm,fig_L,verbose] = create_fig(obj)
verbose      = obj.verbose;
num_workers  = obj.num_workers;
figix        = obj.figix;
dotpm        = obj.dotpm;
run_on_holly = obj.run_on_holly;
distfig      = obj.distfig;

fig_seg = cell(4,1);

if verbose==2 && ~num_workers && ~run_on_holly
    spm_figure('Create','Interactive');
    figure(figix)

    for i=1:size(fig_seg,1)              
        fig_seg{i} = figure(figix + i);
        clf(fig_seg{i})
    end 
else
    for i=1:size(fig_seg,1)              
        close(figure(figix + i));    
    end 
end  

if verbose && dotpm
    fig_L   = figure(figix + 5); clf(fig_L)
    fig_tpm = figure(figix + 6); clf(fig_tpm)    
    if num_workers
        verbose = 0;
    end
else
    fig_tpm = [];
    fig_L   = [];
end

if distfig
    try
        distFig; 
    catch
        warning('distFig not available')
    end
end
drawnow
%==========================================================================

%==========================================================================
function logtpm2tpm(pth_logtpm,dir_write)
if nargin<2, dir_write = './data/results'; end

pth_tpm = 'TPM.nii';
pth_tpm = fullfile(dir_write,pth_tpm);

Nii = nifti(pth_logtpm);
img = Nii.dat(:,:,:,:);
d   = size(img);
dm  = d(1:3);
K   = d(4);

img = reshape(img,[prod(dm) K]);
img = safe_softmax(img);
img = reshape(img,d);

[dir_write,nam,ext] = fileparts(pth_tpm);

vols = cell(K,1);
for k=1:K    
    vols{k}       = fullfile(dir_write,[nam num2str(k) ext]);
    create_nii(vols{k},img(:,:,:,k),Nii.mat,16,nam);
end

matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_tpm;
matlabbatch{1}.spm.util.cat.dtype = 16;
spm_jobman('run',matlabbatch);

delete(fullfile(dir_write,[nam '.mat']));

for k=1:K
    delete(vols{k});
end
%==========================================================================

%==========================================================================
function [tpmnum,tpmden] = smooth_tpm(tpmnum,tpmden,fwhm)
if nargin<3, fwhm = 0.5; end

K = size(tpmnum,4);
for k=1:K
    img = tpmnum(:,:,:,k);            
    img = smooth_img(img,fwhm);
    tpmnum(:,:,:,k) = img;

    img = tpmden(:,:,:,k);            
    img = smooth_img(img,fwhm);
    tpmden(:,:,:,k) = img;
end
%==========================================================================

%==========================================================================
function print_algorithm_progress(status)
fprintf('==============================================\n')  
fprintf(['Algorithm ' status ' ('])
fprintf(datestr(now,'mmmm dd, yyyy HH:MM:SS'))
fprintf(')\n')
fprintf('==============================================\n\n')   
%==========================================================================
