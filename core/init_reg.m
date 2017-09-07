function [pthv,sched,B,a] = init_reg(N,d,nitout,tempdir,Mf)
sched = get_sched(nitout);
sched = [sched(1), sched];

B = affine_basis(12);

f = fullfile(tempdir,'v');
if (exist(f,'dir') == 0)
    mkdir(f);
end

pthv = cell(N,1);
a    = cell(N,1);
for n=1:N
    a{n} = zeros([size(B,3),1]);
    
    v       = zeros([d 3],'single');    
    pthv{n} = fullfile(f,['v' num2str(n) '.nii']);
    create_nii(pthv{n},v,Mf,'float32','Velocity field');
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