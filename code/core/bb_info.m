function [dm,ix_x,ix_y,ix_z] = bb_info(bb)

if numel(bb)==6
    % 3D
    bb = [bb(3) bb(4);bb(1) bb(2);bb(5) bb(6)];
else
    % 2D
    bb = [bb(3) bb(4);bb(1) bb(2);1 1];
end

% Dimensions from bounding box
bb        = floor(bb);
bb(bb==0) = 1;
dm        = diff(bb,1,2)' + 1;

% Indeces from bounding box
ix_x = bb(1,1):bb(1,2);
ix_y = bb(2,1):bb(2,2);
ix_z = bb(3,1):bb(3,2);