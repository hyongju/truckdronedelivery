clear all;close all;clc;
addpath('../lib')

% this script varies the relative efficiency of drones/trucks

global annArborSize 


% annArborSize = 6.6751;
annArborSize = 6.6751;
oneGalJ = 1.3e+8;
oneGalToMile = 10;
oneMileJ = 1.3e+8 / oneGalToMile;
unitDistJ = oneMileJ *annArborSize; 
oneGalPrice = 3.027;
unitPriceUSD = 10/annArborSize * 3.027;

spdCrt = 35 * annArborSize / 60;

nCustomers = 1000;  % NUMBER OF CUSTOMERS' DEMANDS (PICKUP OR DELIVERY)
%nTrucks = [1 2 4 5 8 10 20 25 40];  % default 10
nTrucks = 10;
%maxRad = 0.07;      % MAXIMUM TRAVEL DISTANCE FOR A DRONE BETWEEN A DEMAND LOCATION AND THE HUB(TRUCK) default 0.07
maxRad = 0.4/annArborSize; % default 0.06
relSpeed = 1;     % SPEED OF DRONES / SPEED OF TRUCK
% relEff = 5:5:20;    % RELATIVE EFFICIENCY OF DRONES OVER TRUCK
%relEff = 30;
relEff = [1 5 10 20 30];
ptr = 0;

for u1 = 1:length(relEff)
    for u2 = 1:length(maxRad)
       
        %u2
        [timeCon{u1,u2},feulCon{u1,u2},nDronesCap1{u1,u2},nDronesCap2{u1,u2}] = pathGenAnnArbor(nCustomers,nTrucks,maxRad,relEff(u1),relSpeed,relEff(u1));
        % mTSP (randomized)
        % mTSP (k-mean)
        % mTSP (k-mean) + drones with capacity 1
        % mTSP (k=mean) + drones with capacity 2
        timeCon{u1,u2} = spdCrt * timeCon{u1,u2};
        feulCon{u1,u2} = unitPriceUSD * feulCon{u1,u2};        
    end
end

gmax_x =0;
gmax_y =0;
for i = 1:length(relEff)
    for j = 1:length(maxRad)
        gmax_x = max([max(timeCon{i,j}),gmax_x]);
        gmax_y = max([max(feulCon{i,j}),gmax_y]);
    end
end
gmax_x = 1;gmax_y = 1;
% figure,
% k=1
% for i = 1:length(nTrucks)
%     for j = 1:length(relEff)
%         plot(timeCon{i,j}(k)/gmax_x,feulCon{i,j}(k)/gmax_y,'bs','MarkerSize',10);hold on
%     end
% end
% k=2
% for i = 1:length(nTrucks)
%     for j = 1:length(relEff)
%         plot(timeCon{i,j}(k)/gmax_x,feulCon{i,j}(k)/gmax_y,'ro','MarkerSize',10);hold on
%     end
% end
% k=3
% for i = 1:length(nTrucks)
%     for j = 1:length(relEff)
%         plot(timeCon{i,j}(k)/gmax_x,feulCon{i,j}(k)/gmax_y,'kx','MarkerSize',10);hold on
%     end
% end
% k=4
% for i = 1:length(nTrucks)
%     for j = 1:length(relEff)
%         plot(timeCon{i,j}(k)/gmax_x,feulCon{i,j}(k)/gmax_y,'dm','MarkerSize',10);hold on
%     end
% end

k = 0;
for i = 1:length(relEff)
    for j = 1:length(maxRad)
        k = k+1;
        for l = 1:4
            tC(l,k) = timeCon{i,j}(l)/gmax_x;
            fC(l,k) = feulCon{i,j}(l)/gmax_y;
        end
    end
end
h0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
i=1;
plot(tC(i,:),fC(i,:),'bs','MarkerSize',16,'LineWidth',2);hold on

i=2;
plot(tC(i,:),fC(i,:),'ro','MarkerSize',16,'LineWidth',2);hold on

i=3;
%plot(tC(i,1),fC(i,1),'kh','MarkerSize',16,'LineWidth',3);hold on
plot(tC(i,:),fC(i,:),'kx','MarkerSize',16,'LineWidth',2);hold on
i=4;
plot(tC(i,:),fC(i,:),'dm','MarkerSize',16,'LineWidth',2);hold on
legend('mTSP-RAND','mTSP-K-MEAN','mTSP-K-MEAN-D-C1','mTSP-K-MEAN-D-C2','Location','southeast');


% cset = [];
% for i =1:length(maxRad)
%     cset = [cset max(nDronesCap1{i,1})];
% end
cset2 = relEff;
str = 'travel distance per unit energy UAV to UGV = \{';
for i = 1:length(cset2)
    if i ==length(cset2)
        str = strcat(str,sprintf('%d',cset2(i)),'\}');
    else
        str = strcat(str,sprintf('%d',cset2(i)),',');
    end
end
mxfC = min(50,max(max(fC)));
for j=1:size(tC,2)
    text(tC(4,j),fC(4,j)+0.03*mxfC,num2str(cset2(j)),'FontSize',20,'Color','m');hold on
    text(tC(3,j),fC(3,j)+0.03*mxfC,num2str(cset2(j)),'FontSize',20,'Color','k');hold on
    text(tC(1,j),fC(1,j)+0.03*mxfC,num2str(cset2(j)),'FontSize',20,'Color','b');hold on
    text(tC(2,j),fC(2,j)+0.03*mxfC,num2str(cset2(j)),'FontSize',20,'Color','r');hold on    
end
grid on;
% ang = 0:0.01:pi/2;
% linsp = 0.1:0.1:1.5;
% for j = 1:length(linsp)
%     for i = 1:length(ang)
%         xc = linsp(j)*cos(ang);
%         yc = linsp(j)*sin(ang);
%     end
%     plot(xc,yc,'--','Color',[0.9 0.9 0.9]);hold on;
% end
title(str);
xlabel('average service time (mins)');
ylabel('feul cost (USD)');
% axis([0 1 0 1]);
set(gca,'FontSize',20);
set(findall(h0, 'Type', 'Text'),'FontWeight', 'Normal')
save('allrefEff.mat');