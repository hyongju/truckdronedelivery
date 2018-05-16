function [timeEff,feulEff,maxNumDrone,maxNumDrone2] = pathGenAnnArbor(nCustomers,nTrucks,maxRad,mSpeed,mSpeed2,cntr)
addpath('../lib')
nD = 5;

minNumDrones = 1;
% numLoop =1;
%% INPUT:
% nCustomers: number of customers
% nTrucks: number of trucks
% maxRad: MAXIMUM TRAVEL DISTANCE FOR A DRONE BETWEEN A DEMAND LOCATION AND THE HUB(TRUCK)
% mSpeed: (feul cost of truck / unit distance) / (fuel cost of drone / unit distance) 
% mSpeed2: (speed of drone) / (speed of truck)
%% OUTPUT: 
% timeEff: time efficiency
%       mTSP(random) | mTSP (k-mean) | mTSP (k-mean)&drones (capacity 1) | mTSP (k-mean)&drones (capacity 2)
% feulEff: feul efficiency
%       mTSP(random) | mTSP (k-mean) | mTSP (k-mean)&drones (capacity 1) | mTSP (k-mean)&drones (capacity 2)
% maxNumDrone: max. number of drones used for each truck, when its capacity is 1
% maxNumDrone2: max. number of drones used for each truck, when its capacity is 2

global  min_x min_y max_x max_y

min_x = -83.7996;
max_x = -83.6758;
max_y = 42.3240;
min_y = 42.2227;


% load the border of ann arbor
load('bnd_ann_arbor.mat');

ax_ = b_Poly{2}(:,1);
ay_ = b_Poly{2}(:,2);
% normalize
ax = (ax_ - min(ax_))/(max(ax_)-min(ax_));
ay = (ay_ - min(ay_))/(max(ay_)-min(ay_));


[out1,~]=convhull(ax,ay);
bnd_pnts(:,1) = ax(out1);
bnd_pnts(:,2) = ay(out1);
%plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');

numPts = 0;
tess = convhulln(bnd_pnts);
%testpoints = [ 0 0; 10 10];
% 
% while (numPts < nCustomers)
%     rng('shuffle');
%     rn = rand(1,2);
%     if inhull(rn,bnd_pnts,tess)
%         numPts = numPts + 1
%         pos(numPts,:) = rn;
%     end
%    
% end
load('pos_1000.mat','pos')
%pos = rand(nCustomers,2);



h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);

plot(ax,ay,'k--');hold on;
plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
plot(pos(:,1),pos(:,2),'.');hold on;
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);
set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
xlabel('longitude');ylabel('latitude');
set(gca,'FontSize',16);
axis([0 1 0 1]);axis('equal');

%% mTSP (k-mean)

[clus, cnt] = kmeans(pos,nTrucks);

h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
% 
plot(ax,ay,'k--');hold on;
plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
plot(pos(:,1),pos(:,2),'.');hold on;
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);
set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
xlabel('longitude');ylabel('latitude');
set(gca,'FontSize',16);
axis([0 1 0 1]);axis('equal');

aggCst01 = [];
for i = 1:nTrucks
    cls{i} = pos(find(clus==i),:);
    posc{i} = [0.5 0.5;cls{i};0.5 0.5];
    [~,tour{i}] = tspfunc1(posc{i});
    plot(posc{i}(tour{i},1),posc{i}(tour{i},2));hold on;
%     tLength(i) = calcTourLength(posc{i},tour{i});
    cumSum =0;
    augPth{i} = [];
    for j = 1:length(tour{i})-1
%         tLAug2{i}(1,j) = calcTourLengthNL(pos_cc2{i},tour2{i}(1,1:j+1));
%         posc{i}(tour{i}(j),:)
%         posc{i}(tour{i}(j+1),:)
        [rst_cost,subt{i,j}] = dist2GPS(posc{i}(tour{i}(j),:),posc{i}(tour{i}(j+1),:));
        cumSum = cumSum + rst_cost;
        tLAug{i}(1,j) = cumSum;
        augPth{i} = [augPth{i};subt{i,j}];
    end
    plot(augPth{i}(:,1),augPth{i}(:,2),'-');hold on;
    tLength(i) = cumSum;      
    
    aggCst01 = [aggCst01 tLAug{i}];
%     tLAugT(i) = mean(tLAug{i});
end
totCst01 = sum(aggCst01)/length(aggCst01);


% clf,close
%% mTSP (random)

dvr = nCustomers/nTrucks;

s = RandStream('mlfg6331_64','Seed','shuffle'); % For reproducibility
smpl = datasample(s,1:nCustomers,nCustomers,'Replace',false);

h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);

plot(ax,ay,'k--');hold on;
plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
plot(pos(:,1),pos(:,2),'.');hold on;
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);
set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
xlabel('longitude');ylabel('latitude');
set(gca,'FontSize',16);
axis([0 1 0 1]);axis('equal');

% load('current.mat');

% h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
% 
% plot(ax,ay,'k--');hold on;
% plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
% plot(pos(:,1),pos(:,2),'.');hold on;
% set(gca,'xtick',[0 1]);
% set(gca,'ytick',[0 1]);
% set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
% set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
% xlabel('longitude');ylabel('latitude');
% set(gca,'FontSize',16);
% axis([0 1 0 1]);axis('equal');


for i = 1:nTrucks
    cls1{i} = pos(smpl(dvr*(i-1)+1:dvr*(i-1)+dvr),:);
end
% h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
% plot(pos(:,1),pos(:,2),'.');hold on;
% plot(ax,ay,'k--');hold on;
% plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
% set(gca,'xtick',[0 1]);
% set(gca,'ytick',[0 1]);
% set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
% set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
% set(gca,'FontSize',16);
% xlabel('longitude');ylabel('latitude');
% axis([0 1 0 1]);axis('equal');
aggCst02 = [];
for i = 1:nTrucks
    posc1{i} = [0.5 0.5;cls1{i};0.5 0.5];
    [~,tour1{i}] = tspfunc1(posc1{i});
    plot(posc1{i}(tour1{i},1),posc1{i}(tour1{i},2));hold on;
    augPth1{i} = [];
    cumSum =0;
    for j = 1:length(tour1{i})-1
%         tLAug2{i}(1,j) = calcTourLengthNL(pos_cc2{i},tour2{i}(1,1:j+1));
%         rst_cost = dist2(posc1{i}(tour1{i}(j),:),posc1{i}(tour1{i}(j+1),:));
        [rst_cost,subt1{i,j}] = dist2GPS(posc1{i}(tour1{i}(j),:),posc1{i}(tour1{i}(j+1),:));        
        cumSum = cumSum + rst_cost;
        tLAug1{i}(1,j) = cumSum;
        augPth1{i} = [augPth1{i};subt1{i,j}];
    end
    tLength1(i) = cumSum;    

    
%     for j = 1:length(tour{i})-1
% %         tLAug2{i}(1,j) = calcTourLengthNL(pos_cc2{i},tour2{i}(1,1:j+1));
%         posc{i}(tour{i}(j),:)
%         posc{i}(tour{i}(j+1),:)
%         [rst_cost,subt{i,j}] = dist2GPS(posc{i}(tour{i}(j),:),posc{i}(tour{i}(j+1),:));
%         cumSum = cumSum + rst_cost;
%         tLAug{i}(1,j) = cumSum;
%         augPth{i} = [augPth{i};subt{i,j}];
%     end
    plot(augPth1{i}(:,1),augPth1{i}(:,2),'-');hold on;
%     tLength(i) = cumSum;      
    
    
    
%     
%     %tLength1(i) = calcTourLength(posc1{i},tour1{i});
%     for j = 1:length(tour1{i})-1
%         tLAug1{i}(1,j) = calcTourLengthNL(posc1{i},tour1{i}(1,1:j+1));
%     end
    aggCst02 = [aggCst02 tLAug1{i}];

%     tLAugT1(i) = mean(tLAug1{i});   
end
totCst02 = sum(aggCst02)/length(aggCst02);

% 

%clf,close

% save('current.mat');

%% mTSP (k-mean) + drones with capacity 2

N = minNumDrones* ones(1,nTrucks);
isTrue = zeros(1,nTrucks);


alp = 1.1;
for i =1:nTrucks
%     i
    while(~isTrue(i))
        [cls2{i},cc2{i}] = kmeans(cls{i},N(i));
        rst = [];
        for j = 1:size(cc2{i},1)
            cls2r{i}{j} = cls{i}(find(cls2{i} == j),:);
            if maxRad >= maxDistSet(cc2{i}(j,:),cls2r{i}{j})
                rst = [rst 1];
            else
                rst = [rst 0];
            end
        end
        if all(rst)
            isTrue(i) = 1;
        else
            isTrue(i) = 0;
            N(i) = ceil(N(i) * alp);
        end
    end
end
% h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
% plot(ax,ay,'k--');hold on;
% plot(pos(:,1),pos(:,2),'.');hold on;
% plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
% set(gca,'xtick',[0 1]);
% set(gca,'ytick',[0 1]);
% set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
% set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
% set(gca,'FontSize',16);
% xlabel('longitude');ylabel('latitude');
% axis([0 1 0 1]);axis('equal');

% load('current.mat');

cls2rTmp = cls2r;
for i = 1:nTrucks
    pos_cc2{i} = [ 0.5 0.5;cc2{i};0.5 0.5];
    [~,tour2{i}] = tspfunc1(pos_cc2{i});
    tmlp = pos_cc2{i}(tour2{i},:);
    cc2{i}=tmlp(2:end-1,:);    
%     plot(pos_cc2{i}(tour2{i},1),pos_cc2{i}(tour2{i},2),'LineWidth',3);hold on;
%     tLength2(i) = calcTourLength(pos_cc2{i},tour2{i});
    augPth2{i} = [];
    cumSum =0;
    for j = 1:length(tour2{i})-1
%         tLAug2{i}(1,j) = calcTourLengthNL(pos_cc2{i},tour2{i}(1,1:j+1));
        %rst_cost = dist2(pos_cc2{i}(tour2{i}(j),:),pos_cc2{i}(tour2{i}(j+1),:));
        [rst_cost,subt2{i,j}] = dist2GPS(pos_cc2{i}(tour2{i}(j),:),pos_cc2{i}(tour2{i}(j+1),:));
        cumSum = cumSum + rst_cost;
        tLAug2{i}(1,j) = cumSum;
        augPth2{i} = [augPth2{i};subt2{i,j}];
    end
%     plot(augPth2{i}(:,1),augPth2{i}(:,2),'-');hold on;
    tLength2(i) = cumSum;
    clear idx;
    idx = tour2{i}(1,2:end-1,:);
    for j = 1:size(cc2{i},1)
        
        cls2r{i}{j} = cls2rTmp{i}{idx(j)-1};
    end

%     tLAugT2(i) = mean(tLAug2{i});        
end
% totCst02 = sum(aggCst02)/length(aggCst02);




%clf,close
% save('current.mat');

% mSpeed = 6.1;

h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
plot(ax,ay,'k--');hold on;
plot(pos(:,1),pos(:,2),'.');hold on;
plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);
set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
set(gca,'FontSize',16);
xlabel('longitude');ylabel('latitude');
axis([0 1 0 1]);axis('equal');

maxNumDrone= zeros(1,nTrucks);
accumDcst = zeros(1,nTrucks);


ang = 0:0.01:2*pi;


aggCst = [];

% drone capacity 2
for i = 1:nTrucks
    plot(augPth2{i}(:,1),augPth2{i}(:,2),'-');hold on;
    dLength = 0;
    maxL = 0;
    addafter{i}(1)=0;
    for j = 1:size(cc2{i},1)
        clear tmpvl cst_drone group;
        visited = [];
        ccR{i}{j} = cc2{i}(j,:)+[maxRad*cos(ang)' maxRad*sin(ang)'];
        plot(ccR{i}{j}(:,1),ccR{i}{j}(:,2),'--c');hold on;        
        while(size(visited,2)~=size(cls2r{i}{j},1))
            if length(setdiff(1:size(cls2r{i}{j}),visited)) == 1
                curNode = setdiff(1:size(cls2r{i}{j}),visited);
            else
                curNode = randsample(setdiff(1:size(cls2r{i}{j}),visited),1);
            end
            visited = [visited curNode];
            if size(visited,2)~=size(cls2r{i}{j},1)
                [~,visited] = nearestN(cls2r{i}{j},curNode,visited);
%                 visited = [visited nxtNode];
            end
            %visited = unique(visited);
        end
        %size(visited,2)==size(unique(visited),2)
        l = 0;droneT{i,j} = [];
        for g1 = 1:ceil(size(visited,2)/2)
            l = l+1;
            if (g1==ceil(size(visited,2)/2)) && mod(size(visited,2),2)==1
                group{g1} = cls2r{i}{j}(visited(1,2*(g1-1)+1),:);
                pathD{i,j}{l} = [cc2{i}(j,:);group{g1}(1,:);cc2{i}(j,:)];
                cst_drone(g1) = dist2(cc2{i}(j,:),group{g1}(1,:))*2;
                droneT{i,j} =  [droneT{i,j} dist2(cc2{i}(j,:),group{g1}(1,:))];
            else
                group{g1} = cls2r{i}{j}(visited(1,2*(g1-1)+1:2*g1),:);
                cst_drone(g1) = dist2(cc2{i}(j,:),group{g1}(1,:)) + dist2(cc2{i}(j,:),group{g1}(2,:)) + dist2(group{g1}(1,:),group{g1}(2,:));
                pathD{i,j}{l} = [cc2{i}(j,:);group{g1}(1:2,:);cc2{i}(j,:)];
                droneT{i,j} =  [droneT{i,j} dist2(cc2{i}(j,:),group{g1}(1,:))];
                droneT{i,j} =  [droneT{i,j} dist2(cc2{i}(j,:),group{g1}(1,:))+dist2(group{g1}(1,:),group{g1}(2,:))];
            end
            plot(pathD{i,j}{g1}(:,1),pathD{i,j}{g1}(:,2),'-r');hold on;
             
        end
        numLoop = ceil(l/nD);
        dLengthC{i,j} = sum(cst_drone);
       
        dTimeCst{i,j} = tLAug2{i}(j)+ (numLoop*(ones(size(droneT{i,j}))*mean(droneT{i,j})))/mSpeed2+sum(addafter{i});
        addafter{i}(j) = numLoop*mean(droneT{i,j})/mSpeed2;        
        
        %addafter2{i}(j) = max(droneT2{i,j});
        aggCst = [aggCst dTimeCst{i,j}];
        %dLengthDpT{i,j} = mean(cst_drone) + 
        dLength= dLength + dLengthC{i,j};
        maxL = max(numLoop,maxL);
    end
    accumDcst(i) = accumDcst(i) + dLength;
%         dLength(i) = dLength(i) + sum(sqrt(tmpvl(:,1).^2 + tmpvl(:,2).^2))*2;
    maxNumDrone(i) = ceil(maxL);
end
wdL = accumDcst / mSpeed*2 ;
totCst = sum(aggCst) / length(aggCst);

% axis([0 1 0 1]);axis('equal');

% 
% plot(ax,ay,'k--');hold on;
% plot(pos(:,1),pos(:,2),'.');hold on;
% plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
% set(gca,'xtick',[0 1]);
% set(gca,'ytick',[0 1]);
% set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
% set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
% set(gca,'FontSize',16);
% xlabel('longitude');ylabel('latitude');
% axis([0 1 0 1]);axis('equal');
% clf,close

maxNumDrone2= zeros(1,nTrucks);

%% mTSP (k-mean) + drones with capacity 1

h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
plot(ax,ay,'k--');hold on;
plot(pos(:,1),pos(:,2),'.');hold on;
plot(bnd_pnts(:,1),bnd_pnts(:,2),'k-');
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);
set(gca,'XTickLabel',{sprintf('%.2f',min(ax_)),sprintf('%.2f', max(ax_))});
set(gca,'YTickLabel',{sprintf('%.2f',min(ay_)),sprintf('%.2f', max(ay_))});
set(gca,'FontSize',16);
xlabel('longitude');ylabel('latitude');
axis([0 1 0 1]);axis('equal');


aggCst2 = [];
for i = 1:nTrucks
    plot(augPth2{i}(:,1),augPth2{i}(:,2),'-');hold on;
    dLength2(i) = 0;
    maxL1 = 0;
    addafter2{i}(1)=0;
    % clustering
    for j = 1:size(cc2{i},1)
        clear tmpvl;
        droneT2{i,j} = [];
        tmpvl = repmat(cc2{i}(j,:),size(cls2r{i}{j},1),1) - cls2r{i}{j};
        dLength2(i) = dLength2(i) + sum(sqrt(tmpvl(:,1).^2 + tmpvl(:,2).^2))*2;
        ccR2{i}{j} = cc2{i}(j,:)+[maxRad*cos(ang)' maxRad*sin(ang)'];
        plot(ccR2{i}{j}(:,1),ccR2{i}{j}(:,2),'--c');hold on;             
        for l = 1:size(cls2r{i}{j},1)
            pathD2{i,j}{l} = [cc2{i}(j,:);cls2r{i}{j}(l,:);cc2{i}(j,:)];
            plot(pathD2{i,j}{l}(:,1),pathD2{i,j}{l}(:,2),'-r');hold on;
            droneT2{i,j} =  [droneT2{i,j} dist2(cc2{i}(j,:),cls2r{i}{j}(l,:))];
        end
        numLoop = ceil(l/nD);
        dTimeCst2{i,j} = tLAug2{i}(j)+ (numLoop*(ones(size(droneT2{i,j}))*mean(droneT2{i,j})))/mSpeed2+sum(addafter2{i});
        addafter2{i}(j) = numLoop*mean(droneT2{i,j})/mSpeed2;
        aggCst2 = [aggCst2 dTimeCst2{i,j}];
        maxL1 = max(numLoop,maxL1);
    end
    maxNumDrone2(i) = ceil(maxL1);
    
end
wdL2 = dLength2 / mSpeed;
totCst2 = sum(aggCst2) / length(aggCst2);

% 
%axis([0 1 0 1]);axis('equal');

%clf,close
% h0 = figure('position',[100 100 600 600],'Color',[1 1 1]);
% plot(1:nTrucks,tLength1); hold on;
% plot(1:nTrucks,tLength); hold on;
% plot(1:nTrucks,tLength2+wdL); hold on;
% legend('TSP (random)','TSP (k-mean)','TSP (k-mean)+UAVs');
%clf,close
%out1 = [sum(tLength1),sum(tLength),sum(tLength2),sum(wdL2),sum(wdL)];
timeEff = [totCst02 totCst01 totCst2 totCst];
feulEff = [sum(tLength1),sum(tLength),sum(tLength2)+sum(wdL2),sum(tLength2)+sum(wdL)];

save(sprintf('nTrucks%d.mat',cntr));
