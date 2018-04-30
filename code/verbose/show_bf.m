function show_bf(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       

clims_CT = [0 100];

M   = numel(obj);
cnt = 1;    
for m=1:M      
    modality = obj{m}{1}.modality;
    
    for s=rand_subjs{m} 
        % Get image data
        V   = spm_vol(obj{m}{s}.image(1).fname);
        d0  = V.dim;         
        zix = floor(d0(3)/2) + 1; 
        f   = single(V.private.dat(:,:,zix));       
        
        % Get DCT basis functions
        vx      = spm_misc('vxsize',V.mat);
        fwhm    = obj{m}{s}.segment.biasfwhm(1);        
    
        sd = vx(1)*d0(1)/fwhm; d3(1) = ceil(sd*2);
        sd = vx(2)*d0(2)/fwhm; d3(2) = ceil(sd*2);
        sd = vx(3)*d0(3)/fwhm; d3(3) = ceil(sd*2);
    
        [x0,y0,z0] = ndgrid(1:d0(1),1:d0(2),1);
        B3         = spm_dctmtx(d0(3),d3(3),z0);
        B2         = spm_dctmtx(d0(2),d3(2),y0(1,:)');
        B1         = spm_dctmtx(d0(1),d3(1),x0(:,1));

        % Compute bias field
        T   = obj{m}{s}.segment.Tbias{1};                
        bfz = transf(B1,B2,B3(zix,:),T);
        bfz = exp(bfz);
               
        subplot(8,3,3*(cnt - 1) + 1)
        if strcmp(modality,'CT'), imagesc(f',clims_CT); axis image xy off; title(['f_{' num2str(m), ',' num2str(s), '}']); colormap(gray); colorbar;    
        else                      imagesc(f'); axis image xy off; title(['f_{' num2str(m), ',' num2str(s), '}']); colormap(gray); colorbar;    
        end
        
        subplot(8,3,3*(cnt - 1) + 2)
        imagesc(bfz'); axis image xy off; title(['bf_{' num2str(m), ',' num2str(s), '}']); colormap(gray); colorbar;
                       
        f = bfz.*f;
        
        subplot(8,3,3*(cnt - 1) + 3)
        if strcmp(modality,'CT'), imagesc(f',clims_CT); axis image xy off; title(['bf*f_{' num2str(m), ',' num2str(s), '}']); colormap(gray); colorbar;    
        else                      imagesc(f'); axis image xy off; title(['bf*f_{' num2str(m), ',' num2str(s), '}']); colormap(gray); colorbar;    
        end
        
        cnt = cnt + 1;
    end
end                                  
drawnow
%==========================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%=======================================================================