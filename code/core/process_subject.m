function obj = process_subject(obj,fig) 
try
    rng('default');
    rng(obj.s);
    
    do_template = obj.tot_S>1;

    if obj.segment.est_fwhm && obj.iter==1 && obj.image(1).dim(3)>1
        % Estimate FWHM of image smoothness
        %------------------------------------------------------------------
        N = numel(obj.image);
        fwhm = zeros(1,N);
        for n=1:N
            fwhm(n) = estimate_fwhm(obj.image(n).fname,obj.modality);
        end
        
        obj.fwhm = mean(fwhm);
    end
    
    if (~do_template && ~obj.uniform && obj.maff.do_maff) || ...
       (do_template && obj.maff.maff_done==false && ~obj.uniform && obj.maff.do_maff)
   
        % Affine registration
        %------------------------------------------------------------------
        tpm = spm_load_logpriors(obj.pth_template,obj.segment.wp);

        M                       = obj.image(1).mat;
        c                       = (obj.image(1).dim+1)/2;
        obj.image(1).mat(1:3,4) = -M(1:3,1:3)*c(:);

        [Affine1,ll1]    = spm_maff_new(obj.image(1),8,(0+1)*16,tpm,obj.Affine,obj.maff.affreg);
        Affine1          = Affine1*(obj.image(1).mat/M);
        obj.image(1).mat = M;

        % Run using the origin from the header
        [Affine2,ll2] = spm_maff_new(obj.image(1),8,(0+1)*16,tpm,obj.Affine,obj.maff.affreg);

        % Pick the result with the best fit
        if ll1>ll2, obj.Affine = Affine1; else obj.Affine = Affine2; end

        % Initial affine registration.
        obj.Affine = spm_maff_new(obj.image(1),4,(obj.fwhm+1)*16,tpm,obj.Affine,obj.maff.affreg);            
        obj.Affine = spm_maff_new(obj.image(1),3,obj.fwhm,tpm,obj.Affine,obj.maff.affreg);                            
        clear tpm         
        
        obj.maff.maff_done = true;
        obj.segment.do_def = obj.segment.do_def0;
    elseif do_template && ~obj.uniform && ~obj.maff.do_maff && obj.iter==2
        obj.segment.do_def = obj.segment.do_def0;
    end 

    % Run the new SPM segmentation routine
    %--------------------------------------------------------------  
    tic;    
    obj = spm_segment(obj,fig);
    t1  = toc;    

    % Push responsibilities to template space
    %--------------------------------------------------------------
    t2 = 0;
    if obj.push_resp.do_push_resp
        tic;    
        obj = push_resp(obj);    
        t2  = toc;                            
    end

    % Write results
    %--------------------------------------------------------------
    if obj.write_res.do_write_res      
        write_res(obj);
    end    

    % Verbose
    %--------------------------------------------------------------
    fprintf_obj(obj,t1,t2);   

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
if obj.print_subj_info    
    fprintf('==============================================\n')
    fprintf('%s %s %s %s\n','| modality =',obj.modality,'| fname(1) =',obj.image(1).fname);    
    fprintf('----------------------------------------------\n')
    fprintf('| wp = [%.3f, %s%.3f]\n',obj.segment.wp(1),sprintf('%.3f, ',obj.segment.wp(2:end - 1)),obj.segment.wp(end));
    fprintf('----------------------------------------------\n')
    fprintf('| bf_dc = [%s] | avg_bf_dc = [%s] \n',sprintf('%.3f ', obj.segment.bf_dc),sprintf('%.3f ', obj.segment.avg_bf_dc));
    fprintf('----------------------------------------------\n')    
    for n=1:size(obj.segment.gmm.po.m,1), fprintf('| po.m = [%.3f, %s%.3f]\n',obj.segment.gmm.po.m(n,1),sprintf('%.3f, ',obj.segment.gmm.po.m(n,2:end - 1)),obj.segment.gmm.po.m(n,end)); end
    fprintf('----------------------------------------------\n')
    for n=1:size(obj.segment.gmm.po.m,1), fprintf('| pr.m = [%.3f, %s%.3f]\n',obj.segment.gmm.pr.m(n,1),sprintf('%.3f, ',obj.segment.gmm.pr.m(n,2:end - 1)),obj.segment.gmm.pr.m(n,end)); end
    fprintf('----------------------------------------------\n')
    fprintf('%s %i %s %i %s %i \n','| uniform =',obj.uniform,'| do_bf =',obj.segment.do_bf,'| do_def =',obj.segment.do_def);            
    fprintf('----------------------------------------------\n')
    fprintf('%s %i %s %i %s %.1f %s %.1f %s %.7f \n','| s =',obj.s,'| iter =',obj.iter,'| t1 (s) =',t1,'| t2 (s) =',t2,'| ll =',obj.segment.ll);        
    fprintf('==============================================\n')
end
%==========================================================================