% code for plotting individual trace and its histogram
clear all
close all

% filters variable
filter_rect=8;% s for rho
filter_drift=10; % for leftover drift
frequency=50;%frame rate in Hz

% chose file to analyze
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
addpath(pathname)
name=[pathname, filename];
load(filename);
mkdir([pathname,'figuresrho']);
mkdir([pathname,'figuresrhosquare']);

%make an array of average of x and y (x bar and y bar)
d.bins=filter_rect*frequency+1; % filter 4 s
d.binsdrift=filter_drift*frequency+1;
d.length=e.length;
d.t=e.t;
d.SOURCE=e.SOURCE;

%make an array of average of x and y (x bar and y bar)
for z=1:size(e.x,2)
    d.x_mean(:,z)=smooth(e.x(:,z),d.binsdrift);
    d.y_mean(:,z)=smooth(e.y(:,z),d.binsdrift);
end

%find rho or radius of each bead with respect to the point of attachment
%for drift time window
for z=1:size(e.x,2)
    for k=1:size(e.x,1)
        d.x_cor(k,z)=e.x(k,z)-d.x_mean(k,z);
        d.y_cor(k,z)=e.y(k,z)-d.y_mean(k,z);
        d.x_square(k,z)=d.x_cor(k,z)^2;
        d.y_square(k,z)=d.y_cor(k,z)^2;
        d.rho_square(k,z)=(d.x_cor(k,z)^2)+(d.y_cor(k,z)^2);
        d.rho(k,z)=sqrt(d.rho_square(k,z));
    end
end

% calculating x_square_avg, y_square_avg, rho_square_avg and rho_avg
for z=1:size(d.x_cor,2)
    d.x_square_avg(:,z)=smooth(d.x_square(:,z),d.bins);
    d.y_square_avg(:,z)=smooth(d.y_square(:,z),d.bins);
    d.rho_square_avg(:,z)=smooth(d.rho_square(:,z),d.bins);
    d.rho_avg(:,z)=smooth(d.rho(:,z),d.bins);
end

% defining bins for histograms
i=1;
d.bins_rho_avg(i)=50; %change this according to DNA tether length
while d.bins_rho_avg(i)<425
    d.bins_rho_avg(i+1)=d.bins_rho_avg(i)*1.02;%modify the bin widh manually...here it is 2 percent increase in consecutive bin size
    i=i+1;
end
d.bins_rho_square_avg=d.bins_rho_avg.^2;

%plotting and saving figures for rho
for z=1:size(d.x_cor,2)
    figura=figure(1);
    
    subplot(2,1,1)
    plot(d.t(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),d.rho(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),'*b')
    %         ylim([0 500]) % specify the range manually
    xlabel('time (s)')
    ylabel('<\rho> (nm)')
    title(['bead',num2str(z)]);
    hold on
    plot(d.t(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),d.rho_avg(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),'.r')
    hold off
    
    subplot(2,1,2)
    hhist_rho=histc(d.rho_avg(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),d.bins_rho_avg);
    d.hist_rho_avg(:,z)=hhist_rho;
    hhist_rho=(hhist_rho/(sum(hhist_rho)))./(d.bins_rho_avg.');
    d.pdf_rho_avg(:,z)=hhist_rho;
    bar(d.bins_rho_avg,hhist_rho,'histc')
    %         xlim([100 350])%change limit manually for best view
    xlabel('<\rho> (nm)')
    ylabel('pdf')
    hold off
    
    saveas(figura,[pathname,'figuresrho\','bead ',num2str(z)],'fig')
    clear figura
end

%plotting and saving figures for rhosquare
for z=1:size(d.x_cor,2)
    figura=figure(1);
    
    subplot(2,1,1)
    plot(d.t(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),d.rho_square(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),'*b')
    %         ylim([0 500]) % specify the range manually
    xlabel('time (s)')
    ylabel('<\rho square> (nm^2)')
    title(['bead',num2str(z)]);
    hold on
    plot(d.t(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),d.rho_square_avg(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),'.r')
    hold off
    
    subplot(2,1,2)
    hhist_rho=histc(d.rho_square_avg(round(d.binsdrift/2)+1:d.length(z)-round(d.binsdrift/2),z),d.bins_rho_square_avg);
    d.hist_rho_square_avg(:,z)=hhist_rho;
    hhist_rho=(hhist_rho/(sum(hhist_rho)))./(d.bins_rho_square_avg.');
    d.pdf_rho_square_avg(:,z)=hhist_rho;    
%     hhist_rho=hhist_rho/(sum(hhist_rho));
    bar(d.bins_rho_square_avg,hhist_rho,'histc')
    %         xlim([100 350])%change limit manually for best view
    xlabel('<\rho square> (nm^2)')
    ylabel('pdf')
    hold off
    
    saveas(figura,[pathname,'figuresrhosquare\','bead ',num2str(z)],'fig')
    clear figura
end
close all

%%%%

d.hist_rho_avg_ALL=sum(d.hist_rho_avg,2);
d.hist_rho_square_avg_ALL=sum(d.hist_rho_square_avg,2);

mkdir([pathname,'figuresALL']);

% histogram for all rhoavg
figura=figure(1);
yy=(d.hist_rho_avg_ALL/(sum(d.hist_rho_avg_ALL)))./(d.bins_rho_avg.');
bar(d.bins_rho_avg,yy,'histc')
ylabel('pdf')
xlabel('<\rho> (nm)')
title('ALL beads');

hold off
saveas(figura,[pathname,'figuresALL\','rho_avgALL'],'fig')
clear figura
close all

% histogram for all rhosquareavg
figura=figure(1);
xx=(d.hist_rho_square_avg_ALL/(sum(d.hist_rho_square_avg_ALL)))./(d.bins_rho_square_avg.');
bar(d.bins_rho_square_avg,xx,'histc')
ylabel('pdf')
xlabel('<\rho square> (nm^2)')
title('ALL beads');

hold off
saveas(figura,[pathname,'figuresALL\','rho_square_avgALL'],'fig')
clear figura
close all


%save all the numbers
save('analyzed8s.mat', 'd')% change the name when you change the time window
clear all
close all