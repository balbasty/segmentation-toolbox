function [obj,pars] = init_obj(pars)
M = numel(pars.dat);

%--------------------------------------------------------------------------
dir_subjects = fullfile(pars.dir_output,'subjects');
if exist(dir_subjects,'dir'), rmdir(dir_subjects,'s'); end; mkdir(dir_subjects);   

dir_sum = fullfile(pars.dir_output,'sum');
if exist(dir_sum,'dir'), rmdir(dir_sum,'s'); end; mkdir(dir_sum);   
pth_sgr = fullfile(dir_sum,'sgr.nii');
pth_sH  = fullfile(dir_sum,'sH.nii');

%--------------------------------------------------------------------------
tot_subj = 0;
for m=1:M
    S = numel(pars.dat{m}.V);    
    for s=1:S
        tot_subj = tot_subj + 1; 
    end
end

if ~(tot_subj>1 && pars.do_segment)   
    pars.niter = 1; 
end

%--------------------------------------------------------------------------
fig = cell(1,4);
if pars.dat{1}.segment.verbose
    for i=1:numel(fig), fig{i} = figure(i); clf(figure(i)); end    
end
if pars.niter>1
    if pars.verbose>2, fig{7} = figure(7); clf(figure(7)); end
    if pars.verbose>1, fig{6} = figure(6); clf(figure(6)); end
    if pars.verbose>0, fig{5} = figure(5); clf(figure(5)); end        
end
drawnow

pars.fig = fig;

%--------------------------------------------------------------------------
if M==1, S0 = 8;
else,    S0 = floor(8/M);
end

%--------------------------------------------------------------------------
if pars.do_segment
    V1   = spm_vol(pars.pth_template); 
    d1   = [V1(1).dim numel(V1)];
    d_gr = [d1(1:3),d1(4) - 1];
    d_H  = [d1(1:3) round(((d1(4) - 1)*d1(4))/2)];
end

%--------------------------------------------------------------------------
obj        = cell(1,M);
rand_subjs = cell(1,M);
cnt        = 1;
for m=1:M           
    do_preproc = pars.dat{m}.do_preproc;
    V          = pars.dat{m}.V;
    
    %----------------------------------------------------------------------     
    dir_preproc          = fileparts(V{1}(1).fname);
    dir_preproc          = strsplit(dir_preproc,filesep);
    dir_preproc{end - 1} = [dir_preproc{end - 1} '-preproc'];
    if pars.dat{m}.preproc.rem_neck
        dir_preproc{end - 1} = [dir_preproc{end - 1} '-rn'];
    end
    if pars.dat{m}.preproc.skull_strip
        dir_preproc{end - 1} = [dir_preproc{end - 1} '-ss'];
    end
    dir_preproc = fullfile('/',dir_preproc{2:end - 1});

    if do_preproc   
        if exist(dir_preproc,'dir'), rmdir(dir_preproc,'s'); end; mkdir(dir_preproc);           
    end

    S             = numel(V);
    rand_subjs{m} = randperm(S,min(S0,S));   
    obj{m}        = cell(1,S);
    for s=1:S                
        obj_s = struct;
             
        %------------------------------------------------------------------
        dir_s     = fullfile(dir_subjects,['s-' num2str(s) '_m-' num2str(m)]);
        mkdir(dir_s); 
        obj_s.dir_s = dir_s;
        
        %------------------------------------------------------------------
        N            = numel(V{s});
        obj_s.image  = V{s};
        if isfield(pars.dat{m},'L')
            obj_s.labels = pars.dat{m}.L{s};
        else
            obj_s.labels = [];
        end
        
        %------------------------------------------------------------------
        obj_s.do_segment   = pars.do_segment;
        obj_s.s            = cnt;
        obj_s.pth_template = pars.pth_template;                
        obj_s.tot_S        = tot_subj;
        obj_s.dt           = pars.dt;        
        obj_s.uniform      = pars.uniform;    
        obj_s.dir_template = pars.dir_template;
        obj_s.pth_sgr      = pth_sgr;
        obj_s.pth_sH       = pth_sH;
        obj_s.sum_temp_der = pars.sum_temp_der;
        
        %------------------------------------------------------------------
        obj_s.bf_dc       = zeros(1,N);
        obj_s.avg_bf_dc   = zeros(1,N);
        obj_s.ll_template = 0;
        obj_s.status      = 0;
        obj_s.iter        = 1;
        obj_s.aff_done    = 0;
        
        %------------------------------------------------------------------                              
        dir_s1            = fileparts(V{s}(1).fname);
        dir_s1            = strsplit(dir_s1,filesep);
        dir_s1            = fullfile(dir_preproc,dir_s1{end});                           
        obj_s.dir_preproc = dir_s1;
            
        obj_s.do_preproc                = do_preproc;                                                
        obj_s.preproc.coreg_and_reslice = pars.dat{m}.preproc.coreg_and_reslice;
        obj_s.preproc.do_reslice        = pars.dat{m}.preproc.do_reslice;
        obj_s.preproc.realign2mni       = pars.dat{m}.preproc.realign2mni;
        obj_s.preproc.crop              = pars.dat{m}.preproc.crop;
        obj_s.preproc.crop              = pars.dat{m}.preproc.crop;
        obj_s.preproc.rem_neck          = pars.dat{m}.preproc.rem_neck;
        obj_s.preproc.skull_strip       = pars.dat{m}.preproc.skull_strip;
        
        %------------------------------------------------------------------
        obj_s.bb              = pars.dat{m}.segment.bb;
        obj_s.vox             = pars.dat{m}.segment.vox;
        obj_s.cleanup         = pars.dat{m}.segment.cleanup;
        obj_s.mrf             = pars.dat{m}.segment.mrf;
        obj_s.affreg          = pars.dat{m}.segment.affreg;
        obj_s.wp_lab          = pars.dat{m}.segment.wp_lab;
        obj_s.reg             = pars.dat{m}.segment.reg;
        obj_s.reg0            = obj_s.reg;
        obj_s.fwhm            = pars.dat{m}.segment.fwhm;
        obj_s.samp            = pars.dat{m}.segment.samp;
        obj_s.Affine          = pars.dat{m}.segment.Affine;
        obj_s.biasreg         = pars.dat{m}.segment.biasreg*ones(1,N);
        obj_s.biasfwhm        = pars.dat{m}.segment.biasfwhm*ones(1,N);  
        obj_s.modality        = pars.dat{m}.modality;
        obj_s.do_ml           = pars.dat{m}.segment.do_ml;
        obj_s.do_bf           = pars.dat{m}.segment.do_bf;
        obj_s.do_def          = pars.dat{m}.segment.do_def;
        obj_s.do_wp           = pars.dat{m}.segment.do_wp;
        obj_s.print_seg       = pars.dat{m}.segment.print_seg;
        obj_s.print_ll        = pars.dat{m}.segment.print_ll;        
        obj_s.nsubit          = pars.dat{m}.segment.nsubit;
        obj_s.nitgmm          = pars.dat{m}.segment.nitgmm;
        obj_s.niter           = pars.dat{m}.segment.niter;            
        obj_s.tol1            = pars.dat{m}.segment.tol1;        
        obj_s.kmeans_dist     = pars.dat{m}.segment.kmeans_dist;
        obj_s.init_clust      = pars.dat{m}.segment.init_clust;        
        obj_s.trunc_ct        = pars.dat{m}.segment.trunc_ct;                             
        obj_s.do_missing_data = pars.dat{m}.segment.do_missing_data;
        obj_s.do_write_res    = pars.dat{m}.segment.do_write_res;
        obj_s.do_push_resp    = pars.dat{m}.segment.do_push_resp;
        obj_s.do_old_segment  = pars.dat{m}.segment.do_old_segment;                             
        obj_s.write_tc        = pars.dat{m}.segment.write_tc;
        obj_s.write_bf        = repmat(pars.dat{m}.segment.write_bf,N,1);        
        obj_s.write_df        = pars.dat{m}.segment.write_df;
        
        % Prepare lkp
        %------------------------------------------------------------------
        lkp = pars.dat{m}.segment.lkp;
        
        if isfield(pars.dat{m},'L')
            labels  = pars.dat{m}.L{s};
            labels  = uint8(labels.private.dat(:,:,:));
            mn      = min(labels(labels(:)~=0));
            mx      = max(labels(:));
            
            if isempty(lkp.lab)
                lkp.lab = mn:mx;
                lkp.lab = [lkp.lab zeros(1,pars.K - mx)];
            end
            
%             % Sanity check
%             if max(lkp.lab)>pars.K
%                 error('max(obj_s.lab_lkp)>pars.K')
%             end         
% 
%             if numel(lkp.lab)>numel(mn:mx)
%                 error('numel(obj_s.lab_lkp)>numel(mn:mx)')
%             end
        end
        
        if ~isempty(lkp.rem)
            msk      = ismember(lkp.keep,lkp.rem);
            lkp.keep = lkp.keep(~msk);
        end    
        
        obj_s.lkp = lkp;
        
        %------------------------------------------------------------------
        obj_s.gmm = struct;
        pth_prior = pars.dat{m}.segment.pth_prior;
        if ~isempty(pth_prior) && ~obj_s.do_ml
            tmp          = load(pth_prior,'-mat');
            obj_s.gmm.pr = tmp.pr;
        end
        
        if pars.do_segment
            %------------------------------------------------------------------    
            [~,nam] = fileparts(obj_s.image(1).fname);

            dir_resp = fullfile(dir_s,'resp'); 
            mkdir(dir_resp); 

            pth_resp = cell(1,pars.K);
            for k=1:pars.K
                fname       = fullfile(dir_resp,['resp-', num2str(k), '-' nam, '.nii']);
                pth_resp{k} = fname;  
            end
            obj_s.pth_resp = pth_resp;

            %------------------------------------------------------------------      
            dir_der = fullfile(dir_s,'der');      
            mkdir(dir_der);

            pth_gr = cell(1,d_gr(end));
            for k=1:d_gr(end)
                fname     = fullfile(dir_der,['gr-', num2str(k), '-' nam, '.nii']);
                pth_gr{k} = fname;  
            end
            obj_s.pth_gr = pth_gr;

            pth_H = cell(1,d_H(end));
            for k=1:d_H(end)
                fname    = fullfile(dir_der,['H-', num2str(k), '-' nam, '.nii']);
                pth_H{k} = fname;  
            end
            obj_s.pth_H = pth_H;

            %------------------------------------------------------------------
            dir_vel = fullfile(dir_s,'vel');   
            mkdir(dir_vel);

            obj_s.pth_vel = fullfile(dir_vel,['vel-' nam '.nii']);   

            %------------------------------------------------------------------
            dir_write = fullfile(dir_s,'write');   
            mkdir(dir_write);

            obj_s.dir_write = dir_write;            
        end
        
        %------------------------------------------------------------------
        obj{m}{s} = obj_s;
        
        cnt = cnt + 1;
    end
end

pars.rand_subjs = rand_subjs;
%==========================================================================