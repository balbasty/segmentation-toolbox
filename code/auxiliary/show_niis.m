function show_niis(imdir,show_slices)
if nargin<2, show_slices = true; end

[V,S,N] = get_V(imdir);

fname  = {};
cnt = 1;
for s=1:S
    for n=1:N
        fname{cnt} = V{s}(n).fname;
        cnt        = cnt + 1;    
    end
end

S1    = min(N*S,20);
p     = randperm(N*S,S1);
fname = fname(p);

if show_slices
    figure;
    K1 = floor(sqrt(S1)); K2 = ceil(S1/K1);                                           
    for s=1:S1    
        Nii = nifti(fname{s});
        img = Nii.dat(:,:,:);
        d   = size(img);
        
        subplot(K1,K2,s);
        imagesc(img(:,:,floor(d(3)/2) + 1)'); axis image xy off; colormap(gray);
    end 
    drawnow      
else
    spm_check_registration(char(fname));    
end
