m1 = 1;
s1 = 8;

fprintf('----------------------------------------------\n')
for n=1:size(obj{m1}{s1}.gmm.po.W,1), fprintf('| po.W = [%.7f, %s%.7f]\n',obj{m1}{s1}.gmm.po.W(n,n,1),sprintf('%.7f, ',obj{m1}{s1}.gmm.po.W(n,n,2:end - 1)),obj{m1}{s1}.gmm.po.W(n,n,end)); end
for n=1:size(obj{m1}{s1}.gmm.po.m,1), fprintf('| po.m = [%7f, %s%.7f]\n',obj{m1}{s1}.gmm.po.m(n,1),sprintf('%.7f, ',obj{m1}{s1}.gmm.po.m(n,2:end - 1)),obj{m1}{s1}.gmm.po.m(n,end)); end
fprintf('| po.b = [%.3f, %s%.3f]\n',obj{m1}{s1}.gmm.po.b(1),sprintf('%.3f, ',obj{m1}{s1}.gmm.po.b(2:end - 1)),obj{m1}{s1}.gmm.po.b(end));
fprintf('| po.n = [%.3f, %s%.3f]\n',obj{m1}{s1}.gmm.po.n(1),sprintf('%.3f, ',obj{m1}{s1}.gmm.po.n(2:end - 1)),obj{m1}{s1}.gmm.po.n(end));
fprintf('----------------------------------------------\n')
for n=1:size(obj{m1}{s1}.gmm.pr.W,1), fprintf('| pr.W = [%.7f, %s%.7f]\n',obj{m1}{s1}.gmm.pr.W(n,n,1),sprintf('%.7f, ',obj{m1}{s1}.gmm.pr.W(n,n,2:end - 1)),obj{m1}{s1}.gmm.pr.W(n,n,end)); end
for n=1:size(obj{m1}{s1}.gmm.pr.m,1), fprintf('| pr.m = [%7f, %s%.7f]\n',obj{m1}{s1}.gmm.pr.m(n,1),sprintf('%.7f, ',obj{m1}{s1}.gmm.pr.m(n,2:end - 1)),obj{m1}{s1}.gmm.pr.m(n,end)); end
fprintf('| pr.b = [%.3f, %s%.3f]\n',obj{m1}{s1}.gmm.pr.b(1),sprintf('%.3f, ',obj{m1}{s1}.gmm.pr.b(2:end - 1)),obj{m1}{s1}.gmm.pr.b(end));
fprintf('| pr.n = [%.3f, %s%.3f]\n',obj{m1}{s1}.gmm.pr.n(1),sprintf('%.3f, ',obj{m1}{s1}.gmm.pr.n(2:end - 1)),obj{m1}{s1}.gmm.pr.n(end));
fprintf('----------------------------------------------\n')

obj1              = obj{m1}{s1};
obj1.print_ll     = true;
obj1.print_seg    = true;
obj1.do_push_resp = true;
% obj1.gmm          = rmfield(obj1.gmm,'po');

fig1        = cell(1,4);
for i1=1:numel(fig1), 
    fig1{i1} = figure(i1); clf(figure(i1)); 
end
drawnow

update_subject(obj1,obj1.pth_template,fig1);
