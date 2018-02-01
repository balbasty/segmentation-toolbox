function update_subject(pth_obj,pth_template,new,fig)
obj          = load(pth_obj); 
N            = numel(obj.image);    
write        = obj.write;
do_push_resp = obj.do_push_resp;
bb           = obj.bb;
vx           = obj.vox;
iter         = obj.iter;
tot_S        = obj.tot_S;

build_template = tot_S>1;
tpm            = spm_load_priors8_new(pth_template,build_template);

if ~build_template || iter==2        
    M = obj.image(1).mat;
    c = (obj.image(1).dim+1)/2;
    obj.image(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
    [Affine1,ll1]    = spm_maff8(obj.image(1),8,(0+1)*16,tpm,[],obj.affreg); % Closer to rigid
    Affine1          = Affine1*(obj.image(1).mat/M);
    obj.image(1).mat = M;

    % Run using the origin from the header
    [Affine2,ll2] = spm_maff8(obj.image(1),8,(0+1)*16,tpm,[],obj.affreg); % Closer to rigid

    % Pick the result with the best fit
    if ll1>ll2, obj.Affine  = Affine1; else obj.Affine  = Affine2; end

    % Initial affine registration.
    obj.Affine = spm_maff8(obj.image(1),obj.samp*2,(obj.fwhm+1)*16,tpm, obj.Affine, obj.affreg); % Closer to rigid
    obj.Affine = spm_maff8(obj.image(1),obj.samp*2, obj.fwhm,      tpm, obj.Affine, obj.affreg);
end

reg0 = obj.reg;
if iter>=3
    obj.do_def = true;     
    scal       = 2^max(12 - iter,0);
    obj.reg(3) = obj.reg(3)*scal^2;
end

if new
    % Run the new SPM segmentation algorithm
    tic;    
    obj = spm_preproc8_new(obj,tpm,fig);    
    t   = toc;         
    
    if do_push_resp, 
        ll              = push_resp(obj,tpm,bb,vx); 
        obj.ll_template = ll;
    end
    
    if write, spm_preproc_write8_new(obj,tpm,true(max(obj.lkp),4),true(N,2),true(1,2),obj.mrf,obj.cleanup,obj.bb,obj.vox,obj.dir_write); end
    
    dm_template = size(tpm.dat{1});
    fprintf_obj(obj,t,dm_template);   
else
    % Run the old SPM segmentation algorithm (for comparison)
    tic
    res = spm_preproc8_def(obj);
    fprintf_obj(res,toc);
    
    if 0
        spm_preproc_write8_def(res,true(max(obj.lkp),4),true(N,2),true(1,2),obj.mrf,obj.cleanup,obj.bb,obj.vox,dir_write);
    end
end

obj.reg     = reg0;
obj.uniform = false;
obj.iter    = obj.iter + 1;
save(pth_obj,'-struct','obj')
%==========================================================================

%==========================================================================
function fprintf_obj(obj,t,dm_template)
if true
    fprintf('==============================================\n')
    fprintf('%s %s %s %s\n','| modality =',obj.modality,'| fname(1) =',obj.image(1).fname);    
    fprintf('----------------------------------------------\n')
    fprintf('| wp = [%.3f, %s%.3f]\n',obj.wp(1),sprintf('%.3f, ',obj.wp(2:end - 1)),obj.wp(end));
    fprintf('----------------------------------------------\n')
    fprintf('| bf_dc = [%s] | avg_bf_dc = [%s] \n',sprintf('%.3f ', obj.bf_dc),sprintf('%.3f ', obj.avg_bf_dc));
    fprintf('----------------------------------------------\n')
    if obj.ml
        if obj.new
            for n=1:size(obj.gmm.mn,1), fprintf('| mn = [%.3f, %s%.3f]\n',obj.gmm.mn(n,1),sprintf('%.3f, ',obj.gmm.mn(n,2:end - 1)),obj.gmm.mn(n,end)); end
        else
            for n=1:size(obj.mn,1), fprintf('| mn = [%.3f, %s%.3f]\n',obj.mn(n,1),sprintf('%.3f, ',obj.mn(n,2:end - 1)),obj.mn(n,end)); end
        end
    else    
        for n=1:size(obj.gmm.po.m,1), fprintf('| po.m = [%.3f, %s%.3f]\n',obj.gmm.po.m(n,1),sprintf('%.3f, ',obj.gmm.po.m(n,2:end - 1)),obj.gmm.po.m(n,end)); end
        fprintf('----------------------------------------------\n')
        for n=1:size(obj.gmm.po.m,1), fprintf('| pr.m = [%.3f, %s%.3f]\n',obj.gmm.pr.m(n,1),sprintf('%.3f, ',obj.gmm.pr.m(n,2:end - 1)),obj.gmm.pr.m(n,end)); end
    end
    fprintf('----------------------------------------------\n')
    fprintf('| d0 = [%s] | d1 = [%s] \n',sprintf('%i ', obj.d0),sprintf('%i ', dm_template));        
    fprintf('----------------------------------------------\n')
    fprintf('%s %i %s %i %s %i %s %i %s %i \n','| new =',obj.new,'| ml =',obj.ml,'| uniform =',obj.uniform,'| do_bf =',obj.do_bf,'| do_def =',obj.do_def);            
    fprintf('----------------------------------------------\n')
    fprintf('%s %i %s %i %s %.1f %s %.7f \n','| s =',obj.s,'| iter =',obj.iter,'| time (s) =',t,'| ll =',obj.ll);        
    fprintf('==============================================\n')
end
%==========================================================================