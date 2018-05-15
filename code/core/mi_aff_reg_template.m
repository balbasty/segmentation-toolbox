function obj = mi_aff_reg_template(obj)
if isfield(obj.segment,'wp')
    wp = obj.segment.wp;
else
    V  = spm_vol(obj.pth_template);
    Kb = numel(V);
    wp = ones(1,Kb)/Kb;
end

tpm = spm_load_logpriors(obj.pth_template,wp);

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
%==========================================================================