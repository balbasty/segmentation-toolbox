m = 2;
s = 1;

fprintf('----------------------------------------------\n')
for n=1:size(obj{m}{s}.gmm.po.W,1), fprintf('| po.W = [%.3f, %s%.3f]\n',obj{m}{s}.gmm.po.W(n,n,1),sprintf('%.3f, ',obj{m}{s}.gmm.po.W(n,n,2:end - 1)),obj{m}{s}.gmm.po.W(n,n,end)); end
fprintf('----------------------------------------------\n')
for n=1:size(obj{m}{s}.gmm.po.m,1), fprintf('| po.m = [%.3f, %s%.3f]\n',obj{m}{s}.gmm.po.m(n,1),sprintf('%.3f, ',obj{m}{s}.gmm.po.m(n,2:end - 1)),obj{m}{s}.gmm.po.m(n,end)); end
fprintf('----------------------------------------------\n')
fprintf('| po.b = [%.3f, %s%.3f]\n',obj{m}{s}.gmm.po.b(1),sprintf('%.3f, ',obj{m}{s}.gmm.po.b(2:end - 1)),obj{m}{s}.gmm.po.b(end));
fprintf('----------------------------------------------\n')
fprintf('| po.n = [%.3f, %s%.3f]\n',obj{m}{s}.gmm.po.n(1),sprintf('%.3f, ',obj{m}{s}.gmm.po.n(2:end - 1)),obj{m}{s}.gmm.po.n(end));
fprintf('----------------------------------------------\n')

%%
obj1              = obj{m}{s};
obj1.print_ll     = true;
obj1.print_seg    = true;
obj1.do_push_resp = false;
% obj1.gmm          = rmfield(obj1.gmm,'po');

fig        = cell(1,4);
for i=1:numel(fig), 
    fig{i} = figure(i); clf(figure(i)); 
end
drawnow

update_subject(obj1,pth_template,fig);