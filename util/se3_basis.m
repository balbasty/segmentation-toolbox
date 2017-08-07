function B = se3_basis
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