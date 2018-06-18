function resp = init_resp(obj,lkp,d)

Kb  = max(lkp.part);
mat = obj.image(1).mat;
if ~isfield(obj.segment,'resp')        
    obj.segment.resp = struct;
    if ~isfield(obj.segment.resp,'nii')
        obj.segment.resp.nii = nifti;
        for k=1:Kb
            fname                   = fullfile(obj.dir_seg,['resp-current' num2str(k) '.nii']);    
            obj.segment.resp.nii(k) = spm_misc('create_nii',fname,ones(d,'single')/Kb,mat,spm_type('float32'),'resp-current');
        end
    end   
    
    resp = obj.segment.resp;
else
    resp = obj.segment.resp;
end

d1 = size(obj.segment.resp.nii(1).dat(:,:,:));
if numel(d1)==2, d1 = [d1 1]; end

resp.dat = zeros([d Kb],'single');
if ~isequal(d,d1)
    % Interpolate responsibilities to dimension d
    %---------------------------------------------------------------------            
    s = single(0);
    for k=1:Kb    
        resp.dat(:,:,:,k) = upsample_resp(obj.segment.resp.nii(k).dat(:,:,:),obj); 
        s                 = s + resp.dat(:,:,:,k);
    end
    
    for k=1:Kb    
        resp.dat(:,:,:,k) = resp.dat(:,:,:,k)./s;     
    end
else
    for k=1:Kb    
        resp.dat(:,:,:,k) = single(obj.segment.resp.nii(k).dat(:,:,:));
    end
end
%==========================================================================

%==========================================================================
function img = upsample_resp(img,obj)
d0 = obj.image(1).dim(1:3);
vx = spm_misc('vxsize',obj.image(1).mat);
sk = max([1 1 1],round(obj.segment.samp*[1 1 1]./vx));
    
[x0,y0,z0] = ndgrid(1:d0(1),...
                    1:d0(2),...
                    1:d0(3));

mat0 = eye(4);
mat1 = diag([sk 1]);

T = inv(mat1);
                
x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

img(~isfinite(img))         = 1e-3;
img                         = spm_bsplinc(img,[3 3 3 0 0 0]);
img                         = spm_bsplins(img,x1,y1,z1,[3 3 3 0 0 0]);
img(~isfinite(img) | img<0) = 0;
%==========================================================================