function [timeEff,feulEff,maxNumDrone,maxNumDrone2] = multiRegionPathGen(nCustomers,nTrucks,maxRadSeed,mSpeed,mSpeed2,cntr,bndwh,posD,whloc)
%% this function generates sub-optimal {truck+drones} paths for each region

%%
% INPUT:
% nCustomers: number of customers
% nTrucks: number of trucks
% maxRadSeed: MAXIMUM TRAVEL DISTANCE FOR A DRONE BETWEEN A DEMAND LOCATION AND THE HUB(TRUCK)

% mSpeed: {feul cost of truck / unit distance} / {fuel cost of drone / unit distance} 
% mSpeed2: {speed of drone} / {speed of truck}
% OUTPUT: 
% timeEff: time efficiency
%       mTSP(random) | mTSP (k-mean) | mTSP (k-mean)&drones (capacity 1) | mTSP (k-mean)&drones (capacity 2)
% feulEff: feul efficiency
%       mTSP(random) | mTSP (k-mean) | mTSP (k-mean)&drones (capacity 1) | mTSP (k-mean)&drones (capacity 2)
% maxNumDrone: max. number of drones used for each truck, when its capacity is 1
% maxNumDrone2: max. number of drones used for each truck, when its capacity is 2

nD = 5;                 % number of drones per truck
kmthres = 0.5;          % threshold used for k-mean
minNumDrones = 1;

global  min_x min_y max_x max_y

% border of the region
ax_ = bndwh(:,1);
ay_ = bndwh(:,2);

min_x = min(ax_);
max_x = max(ax_);
min_y = min(ay_);
max_y = max(ay_);

% maximum travel distance of drones (region specific)
maxRad = (max(ay_)-min(ay_))/2.2547 * maxRadSeed;      

% normalize the space
ax = (ax_ - min(ax_))/(max(ax_)-min(ax_))/(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162);
ay = (ay_ - min(ay_))/(max(ay_)-min(ay_));

% normalize warehouse locations
whl(1,1) = (whloc(1) - min(ax_))/(max(ax_)-min(ax_))/(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162);
whl(1,2) = (whloc(2) - min(ay_))/(max(ay_)-min(ay_));



[out1,~]=convhull(ax,ay);
bnd_pnts(:,1) = ax(out1);
bnd_pnts(:,2) = ay(out1);

% normalize demand locations
pos(:,1) = (posD(:,1) - min(ax_))/(max(ax_)-min(ax_))/(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162);
pos(:,2) = (posD(:,2) - min(ay_))/(max(ay_)-min(ay_));

%% mTSP (k-mean-equittable): k-mean clustering where clusters having equal number of points

[clus, ~] = kmeanseq(pos,nTrucks,kmthres);


aggCst01 = [];
for i = 1:nTrucks
    cls{i} = pos(find(clus==i),:);
    posc{i} = [whl;cls{i};whl];
    [~,tour{i}] = tspfunc1(posc{i});
    cumSum =0;
    augPth{i} = [];
    for j = 1:length(tour{i})-1
        [rst_cost,subt{i,j}] = dist2GPS(posc{i}(tour{i}(j),:),posc{i}(tour{i}(j+1),:));
        cumSum = cumSum + rst_cost;
        tLAug{i}(1,j) = cumSum;
        augPth{i} = [augPth{i};subt{i,j}];
    end
    tLength(i) = cumSum;      
    
    aggCst01 = [aggCst01 tLAug{i}];
end
totCst01 = sum(aggCst01)/length(aggCst01);

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


for i = 1:nTrucks
    cls1{i} = pos(smpl(dvr*(i-1)+1:dvr*(i-1)+dvr),:);
end

aggCst02 = [];
for i = 1:nTrucks
    posc1{i} = [whl;cls1{i};whl];
    [~,tour1{i}] = tspfunc1(posc1{i});
    plot(posc1{i}(tour1{i},1),posc1{i}(tour1{i},2));hold on;
    augPth1{i} = [];
    cumSum =0;
    for j = 1:length(tour1{i})-1
        [rst_cost,subt1{i,j}] = dist2GPS(posc1{i}(tour1{i}(j),:),posc1{i}(tour1{i}(j+1),:));        
        cumSum = cumSum + rst_cost;
        tLAug1{i}(1,j) = cumSum;
        augPth1{i} = [augPth1{i};subt1{i,j}];
    end
    tLength1(i) = cumSum;    
    plot(augPth1{i}(:,1),augPth1{i}(:,2),'-');hold on;
    aggCst02 = [aggCst02 tLAug1{i}];

end
totCst02 = sum(aggCst02)/length(aggCst02);


%% mTSP (k-mean) + drones with capacity 2
N = minNumDrones* ones(1,nTrucks);
isTrue = zeros(1,nTrucks);

save(sprintf('newresult31%d.mat',cntr))

alp = 1.1;
for i =1:nTrucks
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
cls2rTmp = cls2r;
for i = 1:nTrucks
    pos_cc2{i} = [whl;cc2{i};whl];
    [~,tour2{i}] = tspfunc1(pos_cc2{i});
    tmlp = pos_cc2{i}(tour2{i},:);
    cc2{i}=tmlp(2:end-1,:);    

    augPth2{i} = [];
    cumSum =0;
    for j = 1:length(tour2{i})-1
        [rst_cost,subt2{i,j}] = dist2GPS(pos_cc2{i}(tour2{i}(j),:),pos_cc2{i}(tour2{i}(j+1),:));
        cumSum = cumSum + rst_cost;
        tLAug2{i}(1,j) = cumSum;
        augPth2{i} = [augPth2{i};subt2{i,j}];
    end

    tLength2(i) = cumSum;
    clear idx;
    idx = tour2{i}(1,2:end-1,:);
    for j = 1:size(cc2{i},1)
        
        cls2r{i}{j} = cls2rTmp{i}{idx(j)-1};
    end
 
end
maxNumDrone= zeros(1,nTrucks);
accumDcst = zeros(1,nTrucks);

ang = 0:0.01:2*pi;

aggCst = [];


for i = 1:nTrucks
    dLength = 0;
    maxL = 0;
    addafter{i}(1)=0;
    for j = 1:size(cc2{i},1)
        clear tmpvl cst_drone group;
        visited = [];
        ccR{i}{j} = cc2{i}(j,:)+[maxRad*cos(ang)' maxRad*sin(ang)'];
        while(size(visited,2)~=size(cls2r{i}{j},1))
            if length(setdiff(1:size(cls2r{i}{j}),visited)) == 1
                curNode = setdiff(1:size(cls2r{i}{j}),visited);
            else
                curNode = randsample(setdiff(1:size(cls2r{i}{j}),visited),1);
            end
            visited = [visited curNode];
            if size(visited,2)~=size(cls2r{i}{j},1)
                [~,visited] = nearestN(cls2r{i}{j},curNode,visited);
            end
        end
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
             
        end
        numLoop = ceil(l/nD);
        dLengthC{i,j} = sum(cst_drone);
       
        dTimeCst{i,j} = tLAug2{i}(j)+ (numLoop*(ones(size(droneT{i,j}))*mean(droneT{i,j})))/mSpeed2+sum(addafter{i});
        addafter{i}(j) = numLoop*mean(droneT{i,j})/mSpeed2;        
        aggCst = [aggCst dTimeCst{i,j}];
        dLength= dLength + dLengthC{i,j};
        maxL = max(numLoop,maxL);
    end
    accumDcst(i) = accumDcst(i) + dLength;
    maxNumDrone(i) = ceil(maxL);
end
wdL = accumDcst / mSpeed*2 ;
totCst = sum(aggCst) / length(aggCst);

maxNumDrone2= zeros(1,nTrucks);

%% mTSP (k-mean) + drones with capacity 1


aggCst2 = [];
for i = 1:nTrucks
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
        for l = 1:size(cls2r{i}{j},1)
            pathD2{i,j}{l} = [cc2{i}(j,:);cls2r{i}{j}(l,:);cc2{i}(j,:)];
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

timeEff = [totCst02 totCst01 totCst2 totCst];
feulEff = [sum(tLength1),sum(tLength),sum(tLength2)+sum(wdL2),sum(tLength2)+sum(wdL)];

save(sprintf('newresult3%d.mat',cntr))
