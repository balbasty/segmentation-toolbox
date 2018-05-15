function write2gif(fig,fname,i)
% Write to gif animation  
frame      = getframe(fig);
im         = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
if i==1
    imwrite(imind,cm,fname,'gif','Loopcount',inf); 
else
    imwrite(imind,cm,fname,'gif','WriteMode','append'); 
end
%==========================================================================