% the beads are selected using clustering and here we plot histogram of mean
% rho square distribution
clear all
close all

[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
load(filename);

clear y
y=[];
z=[];
for j=1:length(d.selected_bead_number)
    i=d.selected_bead_number(j);
    y2=y;
    y1=d.rho_avg(round(d.binsdrift/2)+1:d.length(i)-round(d.binsdrift/2),i);
    z2=z;
    z1=d.rho_square_avg(round(d.binsdrift/2)+1:d.length(i)-round(d.binsdrift/2),i);
    y=cat(1,y1,y2);
    z=cat(1,z1,z2);
end

% making histogram for rho average
[XHIST]=histc(y,d.bins_rho_avg);
d.histSELECTEDrhoavg=XHIST;
d.binsSELECTEDrhoavg=d.bins_rho_avg;

% mean and standard deviation of rho_average of selected beads
d.rhoaverage_selected=y;
d.rhoaverage_selected_mean=mean(y);
d.rhoaverage_selected_std=std(y);

figura=figure(1);  
yy=(XHIST/(sum(XHIST)))./(d.bins_rho_avg.');
bar(d.bins_rho_avg,yy,'histc');
ylabel('pdf')
xlabel('<\rho> (nm)')
title('selected beads');
hold off

saveas(figura,'rho_avg_selected','fig')
clear figura

% making histogram for rho square average
[XHIST]=histc(z,d.bins_rho_square_avg);
d.histSELECTEDrhosquareavg=XHIST;
d.binsSELECTEDrhosquareavg=d.bins_rho_square_avg;

% mean and standard deviation of rho_square_average of selected beads
d.rhosquareavg_selected=z;
d.rhosquareavg_selected_mean=mean(z);
d.rhosquareavg_selected_std=std(z);

figura=figure(1);  
xx=(XHIST/(sum(XHIST)))./(d.bins_rho_square_avg.');
bar(d.bins_rho_square_avg,xx,'histc');
ylabel('pdf')
xlabel('<\rho square> (nm^2)')
title('selected beads');
hold off
% cd(folder_name)
saveas(figura,'rho_square_avg_selected','fig')
clear figura

%save everything
% cd(folder_name)
save('analyzed8s.mat', 'd');
clear all
close all