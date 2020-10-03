% interactively select the rho and rho square traces and plot histograms of
% states
close all
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
load(filename);

% collect the good traces using rho average
clear y
y=[];
m=1;
for i=1:1:size(d.x_cor,2)
    y2=y;
    y1=d.rho_avg(round(d.binsdrift/2)+1:d.length(i)-round(d.binsdrift/2),i);
    
    figura=figure(1);    
    plot(y1,'.b')
    ylim([100 350])
    title(['bead',num2str(i)]);
    xlabel('frames')
    ylabel('<\rho> (nm)')
       
    flag=menu('Will you select this bead for further analysis?','Yes','No');
    if  flag==2
        close all
        clear y1
    elseif flag==1
        y=cat(1,y1,y2);
        f.indices_rhoavg(m)=i;
        m=m+1;
        close all
        clear y1
    end    
end

% making histogram for rho average
[XHIST]=histc(y,d.bins_rho_avg);
f.histSELECTEDrhoavg=XHIST;
f.binsSELECTEDrhoavg=d.bins_rho_avg;
f.rhoaverage=y;

figure(1)
yy=(XHIST/(sum(XHIST)))./(d.bins_rho_avg.');
bar(d.bins_rho_avg,yy,'histc');
ylabel('pdf')
xlabel('<\rho> (nm)')
title('selected beads');
hold off
folder_name = uigetdir;
saveas(figura,[folder_name,'\','rho_avg_selected'],'fig')
clear figura


% collect the good traces using rho square average
clear y
y=[];
m=1;
for i=1:1:size(d.x_cor,2)
    y2=y;
    y1=d.rho_square_avg(round(d.binsdrift/2)+1:d.length(i)-round(d.binsdrift/2),i);
    
    figura=figure(1);    
    plot(y1,'.b')
    ylim([100^2 350^2])
    title(['bead',num2str(i)]);
    xlabel('frames')
    ylabel('<\rho square> (nm^2)')  
        
    flag=menu('Will you select this bead for further analysis?','Yes','No');
    if  flag==2
        close all
        clear y1
    elseif flag==1
        y=cat(1,y1,y2);
        f.indices_rhosquareavg(m)=i;
        m=m+1;
        close all
        clear y1
    end    
end

% making histogram for rho square average
[XHIST]=histc(y,d.bins_rho_square_avg);
f.histSELECTEDrhosquareavg=XHIST;
f.binsSELECTEDrhosquareavg=d.bins_rho_square_avg;
f.rhosquareavg=y;

figure(1)
xx=(XHIST/(sum(XHIST)))./(d.bins_rho_square_avg.');
bar(d.bins_rho_square_avg,xx,'histc');
ylabel('pdf')
xlabel('<\rho square> (nm^2)')
title('selected beads');
hold off
saveas(figura,[folder_name,'\','rho_square_avg_selected'],'fig')
clear figura

%save everything
cd(folder_name)
save('rhoselected.mat', 'f');
clear all
close all