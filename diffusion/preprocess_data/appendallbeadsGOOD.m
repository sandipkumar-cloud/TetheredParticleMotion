%append all the g.mat and h.mat files in a big structure
clear all
glycerol_percent={'0% glycerol','20% glycerol','30% glycerol','50% glycerol','70% glycerol'};
for i=1:5
    flag=menu(glycerol_percent(i),'One');
    [filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
    addpath(pathname)    
    load(filename);
    flag=menu(glycerol_percent(i),'Two');
    [filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
    addpath(pathname)    
    load(filename);
    r.moments(i)=h;
    r.rest(i)=g;
    r.glycerol_percent(i)=glycerol_percent(i);
    clear h g 
end
save('goodbeadsglycerol.mat', 'r');
clear all