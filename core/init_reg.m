function [Pv,sched,B,affpar] = init_reg(N,d,nitout,tempdir,mat)
sched  = get_sched(nitout);  
B      = affine_basis;

f = fullfile(tempdir,'v');
if (exist(f,'dir') == 0)
    mkdir(f);
end

Pv     = cell(N,1);
affpar = cell(N,1);
for n=1:N
    affpar{n} = zeros([size(B,3),1],'single');
    
    v     = zeros([d 3],'single');    
    Pv{n} = fullfile(f,['v' num2str(n) '.nii']);
    create_nii(Pv{n},v,mat,'float32','Velocity field');
end
%==========================================================================

%==========================================================================
function B = affine_basis
B = zeros(4,4,6);
% B = zeros(4,4,12);

B(1,4,1)          = 1;
B(2,4,2)          = 1;
B(3,4,3)          = 1;
B([1,2],[1,2],4)  = [0 1;-1 0];
B([3,1],[3,1],5)  = [0 1;-1 0];
B([2,3],[2,3],6)  = [0 1;-1 0];
B = cat(3,B,diag([1 1 1 0]));
% B(1,1,7)          = 1;
% B(2,2,8)          = 1;
% B(3,3,9)          = 1;
% B([1,2],[1,2],10) = [0 1;1 0];
% B([3,1],[3,1],11) = [0 1;1 0];
% B([2,3],[2,3],12) = [0 1;1 0];
%==========================================================================

%==========================================================================
function sched = get_sched(iter)
shotdef = spm_shoot_defaults;
sched   = shotdef.sched;
if numel(sched) < iter
    sched = [sched, ones(1,iter - numel(sched)) 1];
end
%==========================================================================