function [ax_,ay_]= denormGPS(ax,ay)

global min_x min_y max_x max_y

min_x = -83.7996;
max_x = -83.6758;
max_y = 42.3240;
min_y = 42.2227;


% denormalize
ax_ = ax* (max_x-min_x)+ min_x;
ay_ = ay* (max_y-min_y)+ min_y;
% ay = (ay_ - min(ay_))/(max(ay_)-min(ay_));