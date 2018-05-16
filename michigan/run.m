clear all;close all;clc;
addpath('../lib')
load('../gen_instances/wh_demands.mat');

oneGalJ = 1.3e+8;                   % 1 gallon -- Joule
oneGalToMile = 10;                  % 1 gallon -- mile
oneMileJ = 1.3e+8 / oneGalToMile;   % 1 mile -- Joule
oneGalPrice = 3.027;                % USD per 1 gallon
% design parameters
maxRad = 0.059;                     % MAXIMUM TRAVEL DISTANCE FOR A DRONE BETWEEN A DEMAND LOCATION AND THE HUB(TRUCK) (dimensionless value)
relSpeed = 1;                       % SPEED OF DRONES / SPEED OF TRUCK
relEff = 20;                        % efficiency of drones / trucks
avgSpeedTrucks = 35;                % average truck speed (miles/hrs)
for i = 1:5
    % range of the region in terms of the latitude interval
    citiSize = lldistkm([min_max_xy{i}(1) min_max_xy{i}(3)],[min_max_xy{i}(1) min_max_xy{i}(4)]) / 1.6;
    unitDistJ = oneMileJ *citiSize;     % unit distance -- Joule
    
    % cost to travel unit distance via truck
    unitPriceUSD = citiSize/10 * oneGalPrice;

    % speed conversion (per minute instead of per hour)
    spdCrt = avgSpeedTrucks * citiSize / 60;

    nCustomers = size(posD{i},1);          % size of the demands
    nTrucks = ceil(1*size(posD{i},1)/100); % number of trucks for each region depends on the size of the demand for each region
    
    % this function generates the paths and corresponding costs for 

    % mTSP (randomized)
    % mTSP (k-mean)
    % mTSP (k-mean) + drones with capacity 1
    % mTSP (k-mean) + drones with capacity 2    
    
    [timeCon{i},feulCon{i},nDronesCap1{i},nDronesCap2{i}] = multiRegionPathGen(nCustomers,nTrucks,maxRad,relEff,relSpeed,i,bnd_wh{i},posD{i},C(i,:));
    timeCon{i} = spdCrt * timeCon{i};
    feulCon{i} = unitPriceUSD * feulCon{i};   
    
end