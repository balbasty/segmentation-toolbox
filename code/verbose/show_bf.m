function show_bf(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       
 
M   = numel(obj);
clims_CT = [0 100];
N0       = 1;
for m=1:M
    N = numel(obj{m}{1}.segment.Tbias);
    
    if N>N0
       N0 = N; 
    end
end
 
rows = M*numel(rand_subjs{1});
cols = 3*N0;
 
row = 1;    
for m=1:M      
    modality = obj{m}{1}.modality;
    N        = numel(obj{m}{1}.segment.Tbias);
    
    for s=rand_subjs{m}
        
        for n=1:N
            % Get image data
            V   = spm_vol(obj{m}{s}.image(n).fname);
            d0  = V.dim;         
            zix = floor(d0(3)/2) + 1; 
            f   = single(V.private.dat(:,:,zix));       
 
            % Get DCT basis functions
            vx      = spm_misc('vxsize',V.mat);
            fwhm    = obj{m}{s}.segment.biasfwhm(n);        
 
            sd = vx(1)*d0(1)/fwhm; d3(1) = ceil(sd*2);
            sd = vx(2)*d0(2)/fwhm; d3(2) = ceil(sd*2);
            sd = vx(3)*d0(3)/fwhm; d3(3) = ceil(sd*2);
 
            [x0,y0,z0] = ndgrid(1:d0(1),1:d0(2),zix);
            B3         = spm_dctmtx(d0(3),d3(3),z0);
            B2         = spm_dctmtx(d0(2),d3(2),y0(1,:)');
            B1         = spm_dctmtx(d0(1),d3(1),x0(:,1));
 
            % Compute bias field
            T   = obj{m}{s}.segment.Tbias{n};                
            bfz = transf(B1,B2,B3(zix,:),T);
            bfz = exp(bfz);
 
            col = N0*3*(row - 1);
            
            subplot(rows,cols,col + n)
            if strcmp(modality,'CT'), imagesc(f',clims_CT); axis image xy off; colormap(gray); colorbar;    
            else                      imagesc(f'); axis image xy off; colormap(gray); colorbar;    
            end
 
            subplot(rows,cols,col + N0 + n)
            imagesc(bfz'); axis image xy off;  colormap(gray); colorbar;
 
            f = bfz.*f;
 
            subplot(rows,cols,col + 2*N0 + n)
            if strcmp(modality,'CT'), imagesc(f',clims_CT); axis image xy off; colormap(gray); colorbar;    
            else                      imagesc(f'); axis image xy off; colormap(gray); colorbar;    
            end           
        end
        
        row = row + 1;
    end
end                                  
drawnow
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%==========================================================================