%python_test
%clear all;close all;

function [origM,destM,costM,trajMC] = genRealPath(iOrig, iDest)

kmtomile = 0.621371;
pause(500/1000);
dlmwrite('input_pair.csv',[iOrig iDest], 'delimiter', ',', 'precision', 9);

commandStr = 'python gmaps_io.py';
[status, commandOut] = system(commandStr);
if status==0
 fprintf('squared result is %d\n',str2num(commandOut));
end

M = csvread('result.csv');



flag = zeros(size(M,1),1);
eidx = zeros(size(M,1),1);
for i = 1:size(M,1)
    for j = 5:size(M,2)
        if (j >=6) && (M(i,j) == 0) && (flag(i) ==0)
            eidx(i) = j-1;
            flag(i) = 1;
        end
    end
    if eidx(i) == 0
        eidx(i) = j;
    end
end
eidx = eidx-5;
costM = M(:,5)*kmtomile / 1000;
origM = M(:,1:2);
destM = M(:,3:4);
trajM = M(:,6:end);
for i = 1:size(M,1)
    trajMC{i} = [];
    for j = 1:eidx(i)/2
        trajMC{i} = [trajMC{i};trajM(i,2*(j-1)+1:2*(j-1)+2)];
    end
end
% 
% figure,
% k = 1;
% plot(trajMC{k}(:,2),trajMC{k}(:,1),'-o');hold on;
