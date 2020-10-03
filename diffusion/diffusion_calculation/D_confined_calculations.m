clear all
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);

v.MSD_x_even_avg=mean(v.MSD_x_even,1);
v.MSD_y_even_avg=mean(v.MSD_y_even,1);
v.MSD_x_odd_avg=mean(v.MSD_x_odd,1);
v.MSD_y_odd_avg=mean(v.MSD_y_odd,1);
v.MSD_even_avg=mean(v.MSD_even,1);
v.MSD_odd_avg=mean(v.MSD_odd,1);

clear a b
[a b]=createFit(v.t,v.MSD_x_even_avg);
v.cfit_MSD_x_even_avg=a;
v.struct_MSD_x_even_avg=b;
fig=figure(1);
saveas(fig,'fig_MSD_x_even_avg')
close
v.coeffvalues_MSD_x_even_avg=coeffvalues(a);
v.D_micro_MSD_x_even_avg=v.coeffvalues_MSD_x_even_avg(2)^2*500; % D micro in nanometer square/second
v.L_x_even_avg=v.coeffvalues_MSD_x_even_avg(1);%boundary in nanometers

% clear a b
% [a b]=createFit(v.t,v.MSD_y_even_avg);
% v.cfit_MSD_y_even_avg=a;
% v.struct_MSD_y_even_avg=b;
% fig=figure(1);
% saveas(fig,'fig_MSD_y_even_avg')
% v.coeffvalues_MSD_y_even_avg=coeffvalues(a);
% v.D_micro_MSD_y_even_avg=v.coeffvalues_MSD_y_even_avg(2)^2*500; % D micro in nanometer square/second
% v.L_y_even_avg=v.coeffvalues_MSD_y_even_avg(1);%boundary in nanometers
v.cfit_MSD_x_even= cell(1, size(v.MSD_x_even,1));


for i=1:size(v.MSD_x_even,1)
    clear a b
    [a b]=createFit(v.t,v.MSD_x_even(i,:));
    close
v.cfit_MSD_x_even{i}=a;
v.struct_MSD_x_even(i)=b;

v.coeffvalues_MSD_x_even(i,:)=coeffvalues(a);
v.D_micro_MSD_x_even(i)=v.coeffvalues_MSD_x_even(i,2)^2*500; % D micro in nanometer square/second
v.L_x_even(i)=v.coeffvalues_MSD_x_even(i,1);%boundary in nanometers
end
save('MSD_GOOD.mat','v');
clear all
