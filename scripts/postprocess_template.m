clc; clear;

addpath(genpath('../code'))
addpath('/cherhome/mbrud/dev/distributed-computing')
addpath('/cherhome/mbrud/dev/auxiliary-functions')

%--------------------------------------------------------------------------
%% Paths and dirs
%--------------------------------------------------------------------------

pth_prior    = '/data/mbrud/log-template/res-CHROMIS/prior-CHROMIS-preproc-rn-ss.mat';
pth_template = '/data/mbrud/log-template/res-CHROMIS/template.nii';
dir_out      = '/home/mbrud/Desktop/TEMP/processed-templates/';

%--------------------------------------------------------------------------
%% Make copies of input
%--------------------------------------------------------------------------

% Prior
tmp         = load(pth_prior);
pr0         = tmp.pr;
[~,nam,ext] = fileparts(pth_prior);
pth_nprior  = fullfile(dir_out,['mod-' nam ext]);
copyfile(pth_prior,pth_nprior)

% Template 
[~,nam,ext]   = fileparts(pth_template);
pth_ntemplate = fullfile(dir_out,['mod-' nam ext]);
copyfile(pth_template,pth_ntemplate)

%--------------------------------------------------------------------------
% Post-process...
%--------------------------------------------------------------------------

%% Decrease size
V0 = spm_vol(pth_ntemplate);
bb = [25 140;30 180;40 160];

V1 = cell(1,numel(V0));
for k=1:numel(V0)
    V1{k} = spm_impreproc('subvol',V0(k),bb','sv');        
end

[~,nam,ext]   = fileparts(pth_ntemplate);
pth_ntemplate = fullfile(dir_out,['sv' nam ext]);

%% Modify classes
Nii  = nifti(pth_ntemplate);
img  = Nii.dat(:,:,:,:);
dm   = size(img);
msk  = sum(img(:,:,:,[2 3 5 6 7 8]),4)>log(0.1);
les  = zeros(dm(1:3)) + log(0.2);
les  = msk.*les;
nimg = cat(4,sum(img(:,:,:,[6 7]),4),sum(img(:,:,:,[7 8]),4),img(:,:,:,1),sum(img(:,:,:,[2 3 5]),4),les,img(:,:,:,9));

%% Modify prior
pr = pr0;

pr.b = [mean(pr.b([6 7])) mean(pr.b([7 8]))  pr.b(1) mean(pr.b([2 3 5])) mean(pr.b([2 3 5 6 7 8])) pr.b(9)];
pr.n = [mean(pr.n([6 7])) mean(pr.n([7 8])) pr.n(1) mean(pr.n([2 3 5])) mean(pr.n([2 3 5 6 7 8])) pr.n(9)];
pr.m = [mean(pr.m(:,[6 7]),2) mean(pr.m(:,[7 8]),2) pr.m(:,1) mean(pr.m(:,[2 3 5]),2) mean(pr.m(:,[2 3 5 6 7 8]),2) pr.m(:,9)];
pr.W = cat(3,mean(pr.W(:,:,[6 7]),3),mean(pr.W(:,:,[7 8]),3),pr.W(:,:,1),mean(pr.W(:,:,[2 3 5]),3),mean(pr.W(:,:,[2 3 5 6 7 8]),3),pr.W(:,:,9));

pr.m(:,4) = 100;
pr.m(:,5) = 60;

%% Save new template and prior
[pth,nam,ext] = fileparts(pth_ntemplate);

vols = cell(size(nimg,4),1);
for k=1:size(nimg,4)
    vols{k} = fullfile(pth,[nam num2str(k) ext]);
    create_nii(vols{k},nimg(:,:,:,k),V1{k}.mat,V1{k}.dt,'template');
end

matlabbatch                       = cell(1,1);
matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_ntemplate;    
matlabbatch{1}.spm.util.cat.dtype = 0;
spm_jobman('run',matlabbatch);

delete(fullfile(pth,[nam '.mat']));

for k=1:size(nimg,4)
    delete(vols{k});
end

save(pth_nprior,'pr')