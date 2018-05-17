clear all;close all;clc
% this script simly generate figures showing paths if truck+drones (cap. 2)
% is used

addpath('../lib')

h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);

numRegions = 5;

for pl1 = 1:numRegions
    % load data
    load(sprintf('/results/newresult%d.mat',pl1));
    
    % boundary    
    t1x = (ax - min(ax)) * (max(ax_) - min(ax_))*(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162) + min(ax_);
    t1y = (ay - min(ay)) * (max(ay_) - min(ay_)) + min(ay_);

    plot(t1x,t1y,'k--');hold on;
    
    % demands
    t2x = (pos(:,1) - min(ax)) * (max(ax_) - min(ax_))*(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162) + min(ax_);
    t2y = (pos(:,2) - min(ay)) * (max(ay_) - min(ay_)) + min(ay_);

    plot(t2x,t2y,'k.');hold on;
    
%     t3x = bndwh(:,1);
%     t3y = bndwh(:,2);
%     plot(t3x,t3y,'b--');hold on;

    set(gca,'FontSize',16);
    xlabel('longitude');ylabel('latitude');
    % axis('equal');

    maxNumDrone= zeros(1,nTrucks);
    accumDcst = zeros(1,nTrucks);


    ang = 0:0.01:2*pi;


    aggCst = [];

    % drone capacity 2
    for i = 1:nTrucks
        
        % truck paths
        
        t4x = (augPth2{i}(:,1) - min(ax)) * (max(ax_) - min(ax_))*(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162) + min(ax_);
        t4y = (augPth2{i}(:,2) - min(ay)) * (max(ay_) - min(ay_)) + min(ay_);

        plot(t4x,t4y,'k-');hold on;

        dLength = 0;
        maxL = 0;
        addafter{i}(1)=0;
        for j = 1:size(cc2{i},1)
            clear tmpvl cst_drone group;
            visited = [];
            ccR{i}{j} = cc2{i}(j,:)+[maxRad*cos(ang)' maxRad*sin(ang)'];
    %         plot(ccR{i}{j}(:,1),ccR{i}{j}(:,2),'--m');hold on;        

            t5x = (ccR{i}{j}(:,1) - min(ax)) * (max(ax_) - min(ax_))*(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162) + min(ax_);
            t5y = (ccR{i}{j}(:,2) - min(ay)) * (max(ay_) - min(ay_)) + min(ay_);

%             plot(t5x,t5y,'--c');hold on;        

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
    %             plot(pathD{i,j}{g1}(:,1),pathD{i,j}{g1}(:,2),'-r');hold on;

                t6x = (pathD{i,j}{g1}(:,1) - min(ax)) * (max(ax_) - min(ax_))*(((max(ay_)-min(ay_))/(max(ax_)-min(ax_)))/0.8162) + min(ax_);
                t6y = (pathD{i,j}{g1}(:,2) - min(ay)) * (max(ay_) - min(ay_)) + min(ay_);
                
                % drone paths
                
                plot(t6x,t6y,'-r');hold on;                 

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
end
% overlay google map
plot_google_map('Alpha',0.8,'ShowLabels',0')

