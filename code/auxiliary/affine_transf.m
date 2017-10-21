function varargout = affine_transf(M,varargin)

if numel(varargin)==1
    % Input is one 4D deformation
    varargout{1} = cat(4,M(1,1)*varargin{1}(:,:,:,1) + M(1,2)*varargin{1}(:,:,:,2) + M(1,3)*varargin{1}(:,:,:,3) + M(1,4),...
                         M(2,1)*varargin{1}(:,:,:,1) + M(2,2)*varargin{1}(:,:,:,2) + M(2,3)*varargin{1}(:,:,:,3) + M(2,4),...
                         M(3,1)*varargin{1}(:,:,:,1) + M(3,2)*varargin{1}(:,:,:,2) + M(3,3)*varargin{1}(:,:,:,3) + M(3,4));            
elseif numel(varargin)==3
    % Input are three 3D deformations
    varargout{1} = M(1,1)*varargin{1} + M(1,2)*varargin{2} + M(1,3)*varargin{3} + M(1,4);
    varargout{2} = M(2,1)*varargin{1} + M(2,2)*varargin{2} + M(2,3)*varargin{3} + M(2,4);
    varargout{3} = M(3,1)*varargin{1} + M(3,2)*varargin{2} + M(3,3)*varargin{3} + M(3,4);
else
   error('Unknown input formatting') 
end
%==========================================================================