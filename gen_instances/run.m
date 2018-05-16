clear all;close all;clc
load('data.mat');
% index | population | ----

numW = 5;   % number of warehouses

% 
% for i = 1:size(data,1)
%     cord = [];
%     cord = [cord;data(i,3) data(i,4)];
%     cord = [cord;data(i,5) data(i,6)];
%     cord = [cord;data(i,7) data(i,8)];
%     cord = [cord;data(i,9) data(i,10)];
%     convi{i} = cord;
% end
% figure
% for i = 1:length(convi)
%     K{i} = convhull(convi{i}(:,1),convi{i}(:,2));
%     %plot(convi{i}(K{i},1),-convi{i}(K{i},2),'-');hold on;
%     plot(convi{i}(K{i},1),-convi{i}(K{i},2),'-');hold on;
% end
% 
% 
% 
% numSamp = 5000;
% 
% lb = 1;
% ub = data(1,2);
% lB{1} = lb;
% uB{1} = ub;
% for i = 2:size(data,1)
%     lb = ub + 1;
%     ub = ub + data(i,2);
%     lB{i} = lb;
%     uB{i} = ub;
% end
% for j = 1:numSamp
%     j
%     samp = randsample(1:sum(data(:,2),1),1);
%     whichP = 0;
%     for i = 1:size(data,1)
%         if (samp >= lB{i}) && (samp <= uB{i})
%             whichP = i;
%         end
%     end
%     locsamp(j) = whichP;
% end
% 
% xrange1 = [min(min([data(:,3) data(:,5) data(:,7) data(:,9)])) max(max([data(:,3) data(:,5) data(:,7) data(:,9)]))];
% yrange1 = [min(min([data(:,4) data(:,6) data(:,8) data(:,10)])) max(max([data(:,4) data(:,6) data(:,8) data(:,10)]))];
% 
% rpos = rand(100000,2);
% rposScaled = [rpos(:,1) * xrange1(2) rpos(:,2) * yrange1(2)];
% 
% for i = 1:length(convi)
%     belong2{i} = [];
% end
% 
% for i = 1:size(rposScaled,1)
%     
%     for j = 1:length(convi)
%         polg = [convi{j}(K{j},1) convi{j}(K{j},2)];
%         if inhull(rposScaled(i,:),polg)
%             belong2{j} = [belong2{j} i];
%         end
%     end
% end
% 
% 
% 
% 
% figure,
% for i = 1:length(convi)
%     K{i} = convhull(convi{i}(:,1),convi{i}(:,2));
%     %plot(convi{i}(K{i},1),-convi{i}(K{i},2),'-');hold on;
%     plot(convi{i}(K{i},1),-convi{i}(K{i},2),'-');hold on;
%     for j = 1:length(belong2{i})
%         plot(rposScaled(belong2{i}(j),1),-rposScaled(belong2{i}(j),2),'r.');hold on;
%     end
% end




% figure,
% for i = 1:length(convi)
%     K{i} = convhull(convi{i}(:,1),convi{i}(:,2));
%     %plot(convi{i}(K{i},1),-convi{i}(K{i},2),'-');hold on;
%     plot(convi{i}(K{i},1),-convi{i}(K{i},2),'-');hold on;
% end

% load data to save time...
load('result1.mat');

augSamp = [];
for i = 1:length(locsamp)
    try1 = randsample(belong2{locsamp(i)},1);
    plot(rposScaled(try1,1),-rposScaled(try1,2),'b.'); hold on;
    augSamp = [augSamp;rposScaled(try1,:)];
end

% augSamp = [];
% for i = 1:length(locsamp)
%     augSamp = [augSamp;rposScaled(randsample(belong2{locsamp(i)},1),:)];
% end

% 
%augSamp2 = [augSamp(:,1)-37 -augSamp(:,2)+997];

augSamp2 = [augSamp(:,1) -augSamp(:,2)];

figure, 
plot(augSamp2(:,1),augSamp2(:,2),'.m');

% conversion to GPS locations
lbcorner1 = [-86.824376 41.760299];
rtcorder1 = [-82.506276 43.178177];

lbcorner2 = [37 -997];
rtcorner2 = [794 -659];

augSamp3(:,1) = (augSamp2(:,1) - lbcorner2(1,1))./(rtcorner2(1,1) - lbcorner2(1,1));
augSamp3(:,2) = (augSamp2(:,2) - lbcorner2(1,2))./(rtcorner2(1,2) - lbcorner2(1,2));

augSamp4(:,1) = augSamp3(:,1)*(rtcorder1(1,1) - lbcorner1(1,1)) + lbcorner1(1,1);
augSamp4(:,2) = augSamp3(:,2)*(rtcorder1(1,2) - lbcorner1(1,2)) + lbcorner1(1,2);

h0 = figure('position',[0 0 1300 1200],'Color',[1 1 1]);
%plot(augSamp4(:,1),augSamp4(:,2),'.m'); hold on;
[idx,C] = kmeans(augSamp4, numW);
for j = 1:numW
    posD{j} = augSamp4(find(idx == j),:);
    plot(posD{j}(:,1),posD{j}(:,2),'.'); hold on;
end
plot(C(:,1),C(:,2),'md','MarkerSize',15,'LineWidth',2);hold on;
plot_google_map('Alpha',0.8,'ShowLabels',0')
xlabel('longitude');ylabel('latitude');
set(gca,'FontSize',16);


for j = 1:numW
    h0 = figure('position',[0 0 1300 1200],'Color',[1 1 1]);
    plot(posD{j}(:,1),posD{j}(:,2),'.k'); hold on;
    plot(C(j,1),C(j,2),'md','MarkerSize',15,'LineWidth',2);hold on;
    plot_google_map('Alpha',0.8,'ShowLabels',0')
    xlabel('longitude');ylabel('latitude');
    set(gca,'FontSize',16); 
    pause(1)
end


%% TO JINSUN 
% * numW          number of warehouses
% * C             location of warehouses
% * posD          GPS location of demands associated with warehouses
%                 e.g., posD{i} locations of demands associated with the ith
%                 warehouse


for i = 1:numW
    bnd_wh{i} = posD{i}(convhull(posD{i}(:,1),posD{i}(:,2)),:);
    min_max_xy{i} = [min(min(posD{i}(:,1))) max(max(posD{i}(:,1))) min(min(posD{i}(:,2))) max(max(posD{i}(:,2)))];
end

