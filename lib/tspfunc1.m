function [tour1,tour2] = tspfunc1(pos)
n = size(pos,1);
v_list = 1:n;
e_list = nchoosek(v_list,2);
s = e_list(:,1)';
t = e_list(:,2)';
% wt = rand(1,length(e_list));
for i = 1:length(e_list)
    wt(i) = norm(pos(s(i),:) - pos(t(i),:));
end
for i = 1:length(v_list)
    name{i} = sprintf('%d',i);
end
G = graph(s,t,wt,name);


[T,pred] = minspantree(G,'Type','forest','Root',findnode(G,'1'));
% figure,plot(G,'Layout','layered');


rootedTree = digraph(pred(pred~=0),find(pred~=0),[],G.Nodes.Name);
% figure,plot(rootedTree)
% 
% figure,plot(rootedTree)
% 
% figure,plot(rootedTree)
% 
% figure,plot(rootedTree)
% 
% figure,plot(rootedTree)

e = dfsearch(rootedTree,'1');
for j = 1:length(e)
    tour1(j) = str2num(e{j});
end
disp('size of the input is:');
length(tour1)
[tour2,~] = exchange2(tour1,pos);




