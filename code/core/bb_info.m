function [dm,ix_x,ix_y,ix_z] = bb_info(bb,d)

if numel(bb)==6
    % 3D
    bb = [bb(3) bb(4);bb(1) bb(2);bb(5) bb(6)];
else
    % 2D
    bb = [bb(3) bb(4);bb(1) bb(2);1 1];
end

bb      = bb';
bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
bb(2,:) = min(bb(2,:),d);

dm = diff(bb) + 1;

% Indeces from bounding box
ix_x = bb(1,1):bb(2,1);
ix_y = bb(1,2):bb(2,2);
ix_z = bb(1,3):bb(2,3);