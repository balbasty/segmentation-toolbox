function [obj,munum1,muden1] = update_model_parameters(obj,logtpm,iter)

if obj.doaff && ((obj.use_tpm && iter==1) || (~obj.use_tpm && iter==3))
    % Affinely register template to 1st image
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,logtpm,obj.Affine,'mni');            
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp, obj.fwhm,      logtpm,obj.Affine,'mni');                                    
end                 

if obj.dodef0 && obj.use_tpm && iter==1
    % If a template is used, start estimating deformations from 1st
    % iteration
    obj.dodef = 1;
end

if ~obj.use_tpm && iter==1 && obj.dotpm   
    % 1st iteration
    obj.lkp   = 1:max(obj.lkp); % One Gaussian per tissue
    obj.dodef = 0;              % No deformation update for 1st iteration
elseif ~obj.use_tpm && iter==2 && obj.dotpm                      
    % 2nd iteration
    obj.lkp   = repelem(1:max(obj.lkp),obj.nlkp); % Use nlkp Gaussians per tissue

    % Re-estimate cluster parameters based on template constructed after 1st iteration
    if isfield(obj,'mg'), obj = rmfield(obj,'mg'); end
    if isfield(obj,'wp'), obj = rmfield(obj,'wp'); end
    if isfield(obj,'mn'), obj = rmfield(obj,'mn'); end
    if isfield(obj,'vr'), obj = rmfield(obj,'vr'); end
    if isfield(obj,'po'), obj = rmfield(obj,'po'); end
    if isfield(obj,'pr'), obj = rmfield(obj,'pr'); end        
elseif ~obj.use_tpm && iter==4 && obj.dotpm       
    % 4th iteration
    obj.dodef      = obj.dodef0; % Start estimating deformations
    obj.niter_stop = 1;          % Allow algorithm to stop earlier, if converged
end         

if isfield(obj,'pthTwarp')
    % To save memory, load deformation from file
    load(obj.pthTwarp,'Twarp')
    obj.Twarp = Twarp;
end

% Run segmentation algorithm
try
    [obj,munum1,muden1] = spm_preprocx(obj,logtpm);
catch
    warning('spm_preprocx')

    obj.ll = 0;
    obj.nm = 0;
    munum1 = single(0);
    muden1 = single(0);    
end

if isfield(obj,'pthTwarp')
    % To save memory, delete deformation from result structure
    Twarp = obj.Twarp;
    save(obj.pthTwarp,'Twarp')
    obj   = rmfield(obj,'Twarp');
end
