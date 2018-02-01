function run_segment

% -Push to GitHub
% -Work on CT
% -Add missing (if artefact, try not doing circular bsplines)
% -Introduce aux
% -Parameterise in log-space

addpath(genpath('./spm-default'))
addpath(genpath('./spm-new'))

%--------------------------------------------------------------------------
pars.print_ll = false;
pars.show_fig = false;

pars.ml     = false;
pars.do_bf  = true;
pars.do_def = true;
pars.write  = false;

pth_template    = '';
% pth_template    = fullfile(spm('dir'),'tpm','TPM.nii');  

pars.dir_output = '/home/mbrud/Desktop/build-template-output/';

K = 6;
S = 16; 

%-------------------------------------------------------------------------- 
im = {};
im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/IXI-noneck/',S,'MRI',0,3,''}; 
% im{end + 1} = {fullfile(get_pth_dropbox,'/PhD/Data/IXI'),S,'MRI',0,4,''}; 
% im{end + 1} = {fullfile(get_pth_dropbox,'/PhD/Data/CHROMIS'),S,'CT',0,4,''}; 
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-CHROMIS-noneck/',S,'CT',0,3,''};  

%--------------------------------------------------------------------------
M               = numel(im);
V               = cell(1,M);
for m=1:M, V{m} = get_V(im{m}); end

%--------------------------------------------------------------------------
[pth_template,uniform] = init_uniform_template(pth_template,V,K,pars.dir_output); 

%------------------------------------------------------
[pth_obj,niter,fig,rand_subjs,pth_H,pth_gr] = init_obj(V,im,pth_template,uniform,pars);
      
%--------------------------------------------------------------------------
print_algorithm_progress('started')

L = -Inf;
for iter=1:niter        

    %----------------------------------------------------------------------    
    reset_template_derivatives(pth_obj,pth_template,pth_H,pth_gr);
    
    % Compute average bias field DC component from pth_obj
    %----------------------------------------------------------------------    
    get_avg_bf_dc(pth_obj);       
    
    % Segment a bunch of subjects
    %----------------------------------------------------------------------
    for m=1:numel(pth_obj)
        S         = numel(pth_obj{m});   
        pth_obj_m = pth_obj{m};
        if S==1
            for s=1:S
                update_subject(pth_obj_m{s},pth_template,true,fig);
            end    
        else
            parfor s=1:S
                update_subject(pth_obj_m{s},pth_template,true,fig);
            end
        end
    end
    
    if niter>1
        % Update template
        %----------------------------------------------------------------------
        L = update_template(L,pth_template,pth_obj,pth_H,pth_gr,iter);                
        
        % Show some output
        show_resp(fig{5},pth_obj,rand_subjs);        
        show_template(fig{6},pth_template,pth_obj);        
        plot_ll(L,fig{7})
    end
    
    if niter>1
        % Update Gaussian-Wishart hyper-parameters
        %----------------------------------------------------------------------
        update_intensity_prior(pth_obj);
    end    
end

print_algorithm_progress('finished')
%==========================================================================

%==========================================================================
function reset_template_derivatives(pth_obj,pth_template,pth_H,pth_gr)
for m=1:numel(pth_obj)
    S = numel(pth_obj{m}); 
    for s=1:S
        obj             = load(pth_obj{m}{s}); 
        obj.ll_template = 0;
        save(pth_obj{m}{s},'-struct','obj') 
    end    
end

V = spm_vol(pth_template);
d = [V(1).dim numel(V)];
create_nii(pth_gr,zeros([d(1:3),d(4)-1],'single'),eye(4),'float32','gr');    
create_nii(pth_H,zeros([d(1:3) round(((d(4)-1)*d(4))/2)],'single'),eye(4),'float32','H');   
%==========================================================================

%==========================================================================
function get_avg_bf_dc(pth_obj)
M = numel(pth_obj);
for m=1:M
    S         = numel(pth_obj{m});
    sum_bf_dc = 0;
    for s=1:S
        obj       = load(pth_obj{m}{s},'-mat','bf_dc');  
        sum_bf_dc = sum_bf_dc + obj.bf_dc;
    end
    clear obj
    
    avg_bf_dc = sum_bf_dc/S;
    
    for s=1:S
        obj           = load(pth_obj{m}{s});    
        obj.avg_bf_dc = avg_bf_dc;
        save(pth_obj{m}{s},'-struct','obj')   
    end
end
%==========================================================================

%==========================================================================
function show_template(fig,pth_template,pth_obj)
obj = load(pth_obj{1}{1},'gmm');
pr  = obj.gmm.pr;

set(0,'CurrentFigure',fig);     

Nii = nifti(pth_template);
b   = Nii.dat(:,:,:,:);
d   = size(b);
K   = d(4);                                  
for k=1:K         
    pr_m1 = num2str(round(pr.m(1,k),2));
    
    subplot(3,K,k);
    slice = b(:,:,floor(d(3)/2) + 1,k);
    imagesc(slice'); axis image xy off; title(['k=' num2str(k)]); colormap(pink);               
   
    subplot(3,K,K + k);
    slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
    imagesc(slice); axis image xy off; title(['m=' pr_m1]); colormap(pink);   

    subplot(3,K,2*K + k);
    slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
    imagesc(slice'); axis image xy off; colormap(pink);   
end 
drawnow
%==========================================================================

%==========================================================================
function show_resp(fig,pth_obj,rand_subjs)
M = numel(pth_obj);
set(0,'CurrentFigure',fig);       

cnt = 1;    
for m=1:M                
    for s=rand_subjs{m}
        obj = load(pth_obj{m}{s},'-mat');    
        K   = numel(obj.pth_resp);
        
        for k=1:K    
            Nii = nifti(obj.pth_resp{k});
            img = Nii.dat(:,:,:);
            dm  = size(img);
            zix = floor(dm(3)/2) + 1;
            img = img(:,:,zix);
        
            subplot(M*numel(rand_subjs{1}),K,cnt);
            imagesc(img'); axis image xy off; colormap(pink);
            title(['q_{' num2str(m), ',' num2str(s) ',' num2str(k) '}']);
            cnt = cnt + 1;
        end 
    end
end                                  
drawnow
%==========================================================================

%==========================================================================
function print_algorithm_progress(status)
fprintf('==============================================\n')  
fprintf(['| Algorithm ' status ' ('])
fprintf(datestr(now,'mmmm dd, yyyy HH:MM:SS'))
fprintf(')\n')
fprintf('==============================================\n\n')   
%==========================================================================

%==========================================================================
function plot_ll(L,fig)
set(0,'CurrentFigure',fig);                
plot(0:numel(L(3:end)) - 1,L(3:end),'b-','LineWidth',1);   hold on            
plot(0:numel(L(3:end)) - 1,L(3:end),'b.','markersize',10); hold off  
title('ll')
%==========================================================================