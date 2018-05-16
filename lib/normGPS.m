function [ax,ay]= normGPS(ax_,ay_)

global min_x min_y max_x max_y

min_x = -83.7996;
max_x = -83.6758;
max_y = 42.3240;
min_y = 42.2227;


% min_x = min(ax_);
% min_y = min(ay_);
% max_x = max(ax_);
% max_y = max(ay_);

% normalize
ax = (ax_ - min_x)/(max_x-min_x);
ay = (ay_ - min_y)/(max_y-min_y);

