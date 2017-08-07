function [pthv,sched,B,affpar] = init_reg(N,d,nitout,tempdir,mat)
sched = get_sched(nitout);  
B     = se3_basis;

f = fullfile(tempdir,'v');
if (exist(f,'dir') == 0)
    mkdir(f);
end

pthv     = cell(N,1);
affpar = cell(N,1);
for n=1:N
    affpar{n} = zeros([size(B,3),1],'single');
    
    v     = zeros([d 3],'single');    
    pthv{n} = fullfile(f,['v' num2str(n) '.nii']);
    create_nii(pthv{n},v,mat,'float32','Velocity field');
end
%==========================================================================

%==========================================================================
function sched = get_sched(iter)
shotdef = spm_shoot_defaults;
sched   = shotdef.sched;
if numel(sched) < iter
    sched = [sched, ones(1,iter - numel(sched)) 1];
end
%==========================================================================