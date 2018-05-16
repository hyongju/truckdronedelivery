function [out1,out2] = dist2GPS(p1,p2)
global annArborSize

% clear all;close all;

% annArborSize = 6.8

% p1= [34.042715	-118.266192]
% p2 = [42.273	-83.738]

% n = size(pos,1);
% v_list = 1:n;
% e_list = nchoosek(v_list,2);
% s = e_list(:,1)';
% t = e_list(:,2)';
% % wt = rand(1,length(e_list));
% length(e_list)
% for i = 1:length(e_list)
%     i
% %     wt(i) = norm(pos(s(i),:) - pos(t(i),:));

[iOrigx,iOrigy] = denormGPS(p1(1),p1(2));
[iDestx,iDesty] = denormGPS(p2(1),p2(2));

iOrig = [iOrigy,iOrigx];
iDest = [iDesty,iDestx];
[~,~,costM,trajM] = genRealPath(iOrig, iDest);
out1 = costM/6.8;
out2Tmp = [iOrig;trajM{1};iDest];

for i = 1:size(out2Tmp,1)
    [out2(i,1),out2(i,2)]= normGPS(out2Tmp(i,2),out2Tmp(i,1));
    
end

%out2
% 
% for i = 1:length(v_list)
%     name{i} = sprintf('%d',i);
% end

% cPath = [];
% wtt = 0;
% for i = 1:length(tour2)-1
%     a1 = find(e_list(:,1) == tour2(i))';
%     a2 = find(e_list(:,2) == tour2(i+1))';
%     edj = intersect(a1,a2);
%     if isempty(edj)
%         a1 = find(e_list(:,1) == tour2(i+1))';
%         a2 = find(e_list(:,2) == tour2(i))';        
%         edj = intersect(a1,a2);
%         cPath = [cPath;flip(trajMC{edj}{1},1)];
%     else
%         cPath = [cPath;trajMC{edj}{1}];
%     end
% %     cPath = [cPath;trajMC{edj}{1}];
%     wtt = wtt + wt(edj) ;
% end


