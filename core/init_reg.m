function [pthv,sched,B,a,R,prm0,prm,int_args,vconv] = init_reg(N,d,nitmain,tempdir,Mf,abasis,alam,vx,vlam)

B = affine_basis(abasis);

R(1:3,1:3)     = alam*eye(3); % translation
R(4:6,4:6)     = alam*eye(3); % rotation
R(7:9,7:9)     = alam*eye(3); % scaling
R(10:12,10:12) = alam*eye(3); % skew    

sched = get_sched(nitmain);
prm0  = [vx vlam*sched(1)];
prm   = prm0;

int_args = 1;

vconv = ones(1,N,'logical'); % Keeps track of convergence of velocity fields, so that the correct regularisation is used

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
function sched = get_sched(nitmain)
shotdef = spm_shoot_defaults;
sched   = shotdef.sched;
sched   = [sched(2) sched(2:end)];
if numel(sched) < nitmain
    sched = [sched, ones(1,nitmain - numel(sched)) 1];
end
%==========================================================================