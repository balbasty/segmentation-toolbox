function obj = process_subject(obj,fig) 
try
    rng('default');
    rng(obj.s);
    
    do_template = obj.tot_S>1;

    if (~do_template && ~obj.uniform && obj.do_aff) || ...
       (do_template && obj.aff_done==false && ~obj.uniform && obj.do_aff)
        % Affine registration
        %--------------------------------------------------------------
        tpm = spm_load_logpriors(obj.pth_template);

        M                       = obj.image(1).mat;
        c                       = (obj.image(1).dim+1)/2;
        obj.image(1).mat(1:3,4) = -M(1:3,1:3)*c(:);

        if obj.image(1).dim(3)>1
            % Input image(s) is 3D
            [Affine1,ll1]    = spm_maff_new(obj.image(1),8,(0+1)*16,tpm,obj.Affine,obj.affreg);
            Affine1          = Affine1*(obj.image(1).mat/M);
            obj.image(1).mat = M;

            % Run using the origin from the header
            [Affine2,ll2] = spm_maff_new(obj.image(1),8,(0+1)*16,tpm,obj.Affine,obj.affreg);

            % Pick the result with the best fit
            if ll1>ll2, obj.Affine = Affine1; else obj.Affine = Affine2; end

            % Initial affine registration.
            obj.Affine = spm_maff_new(obj.image(1),4,(obj.fwhm+1)*16,tpm,obj.Affine,obj.affreg);            
            obj.Affine = spm_maff_new(obj.image(1),3,obj.fwhm,tpm,obj.Affine,obj.affreg);                    
        else
            % Input image(s) is 2D
            obj.Affine = spm_maff_new(obj.image(1),2,obj.fwhm,tpm,obj.Affine,'rigid');  
        end
        clear tpm 

        obj.aff_done = true;                        
        obj.do_def   = true; % activate non-linear registration
    elseif do_template && ~obj.uniform && ~obj.do_aff && obj.iter==2
        obj.do_def   = true; % activate non-linear registration
    end 

    if ~obj.do_old_segment
        % Run the new SPM segmentation routine
        %--------------------------------------------------------------  
        tic;    
        obj = spm_segment(obj,fig);
        t1  = toc;    

        % Push responsibilities to template space
        %--------------------------------------------------------------
        t2 = 0;
        if obj.do_push_resp
            tic;    
            obj = push_resp(obj);    
            t2  = toc;                            
        end

        % Write results
        %--------------------------------------------------------------
        if obj.do_write_res      
            write_res(obj);
        end    

        % Verbose
        %--------------------------------------------------------------
        fprintf_obj(obj,t1,t2);   

    else
        % Run the old SPM segmentation algorithm (for comparison)
        %--------------------------------------------------------------
        res = spm_preproc8_def(obj);        

        spm_preproc_write8_def(res,obj.write_tc,obj.write_bf,obj.write_df,obj.mrf,obj.cleanup,obj.bb,obj.vox,obj.dir_write);
    end

    obj.status = 0; % success
catch ME            
    fprintf(['Error for image: ' obj.image(1).fname '\n'])
    for i=1:numel(ME.stack)
        disp([ME.stack(i).name ', line ' num2str(ME.stack(i).line)]);
    end
    disp(ME.message)    

    obj.status = 1; % fail
end
%==========================================================================

%==========================================================================
function fprintf_obj(obj,t1,t2)
if obj.print_seg
    V  = spm_vol(obj.pth_template);
    d1 = V(1).dim;

    fprintf('==============================================\n')
    fprintf('%s %s %s %s\n','| modality =',obj.modality,'| fname(1) =',obj.image(1).fname);    
    fprintf('----------------------------------------------\n')
    fprintf('| wp = [%.3f, %s%.3f]\n',obj.wp(1),sprintf('%.3f, ',obj.wp(2:end - 1)),obj.wp(end));
    fprintf('----------------------------------------------\n')
    fprintf('| bf_dc = [%s] | avg_bf_dc = [%s] \n',sprintf('%.3f ', obj.bf_dc),sprintf('%.3f ', obj.avg_bf_dc));
    fprintf('----------------------------------------------\n')
    if obj.do_ml
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
    fprintf('| d0 = [%s] | d1 = [%s] \n',sprintf('%i ', obj.d0),sprintf('%i ', d1));        
    fprintf('----------------------------------------------\n')
    fprintf('%s %i %s %i %s %i %s %i %s %i \n','| new =',obj.new,'| ml =',obj.do_ml,'| uniform =',obj.uniform,'| do_bf =',obj.do_bf,'| do_def =',obj.do_def);            
    fprintf('----------------------------------------------\n')
    fprintf('%s %i %s %i %s %.1f %s %.1f %s %.7f \n','| s =',obj.s,'| iter =',obj.iter,'| t1 (s) =',t1,'| t2 (s) =',t2,'| ll =',obj.ll);        
    fprintf('==============================================\n')
end
%==========================================================================