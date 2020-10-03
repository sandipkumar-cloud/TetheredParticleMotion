%TPM data analysis CONTROL $interlaced files$ $50Hz$ $filter rectangular 4 s$
%original file by Carlo Manzo for version 2006b
%modifications by David Dunlap for version 2011a
%modification by Sandip Kumar for variable time filters: version 2011a
%modification by Sandip Kumar for  a different calculation of rho: version
%2011a
%correction of the pixel size
%incorporation of better symmetry test.. covariance matrix (lin  han et al.)
%http://en.wikipedia.org/wiki/Covariance_matrix
%above link explains why we take the ratio of the square root of the eigen
%values(square root of the eigen values represent the lengths of the major and minor axes)
%plot radial histogram

clear all

close all
% filters variable
filter_drift=40;% s
filter_rect=[0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12,10.24,20.48,40.96,81.92];
bins=(filter_rect/0.02)+1;


%%%%%%%%%%%%
% chose file to analyze
ext='.txt';
[filename, pathname] = uigetfile({'*.txt';'*.*'},'File Selector');

%remove figures folder and make a new one
if exist([pathname,'figures'],'dir')
    rmdir([pathname,'figures'],'s');
end
mkdir([pathname,'figures']);

%remove all old .mat file
delete([pathname,'*.mat']);

%for 100X, pixel size is 65nm and for 63X, pixel size is 103.17nm and.
pix_size_1=(6500/63);
pix_size_2=(6250/63);


% %%%%%%%read and show the data

a.data=dlmread([pathname,filename],'',3,0);%reads the text file into structure a
%a skips the first two lines of data as there is a break after this. so we
%start reading from data line 3


%odd and even lines are separated and the data is stored in b. the
%separation is required because of the difference in the odd and even
%sensors in the camera
b.time_odd=(a.data(1:2:length(a.data(:,1)),1)-a.data(1,1))/1000;
b.data_odd=(a.data(1:2:length(a.data(:,1)),2:size(a.data,2)));
b.time_even=(a.data(2:2:length(a.data(:,1)),1)-a.data(1,1))/1000;
b.data_even=(a.data(2:2:length(a.data(:,1)),2:size(a.data,2)));
b.beadsnumber=size(b.data_odd,2)/3;
clear a

%coordinates are converted to nm or length scale
for i=1:b.beadsnumber
    b.beadslabel(i)=b.data_odd(1,i*3-2);
    b.beadsx_odd(:,i)=b.data_odd(:,i*3-1)*pix_size_1;
    b.beadsy_odd(:,i)=b.data_odd(:,i*3)*pix_size_2;
    
    b.beadsx_even(:,i)=b.data_even(:,i*3-1)*pix_size_1;
    b.beadsy_even(:,i)=b.data_even(:,i*3)*pix_size_2;
end


%the coordinates of all the beads are displayed. the even and odd frames
%are shown separately

%odd frames
figure('position',[10 550 560 420]);
figure(1)
plot(b.beadsx_odd,b.beadsy_odd)
axis square equal
title('Odd frames, beads image')
xlabel('x position (nm)')
ylabel('y position (nm)')
if b.beadsnumber==1
    legend('bead 1')
elseif b.beadsnumber==2
    legend('bead 1', 'bead 2')
elseif b.beadsnumber==3
    legend('bead 1', 'bead 2', 'bead 3')
elseif b.beadsnumber==4
    legend('bead 1', 'bead 2', 'bead 3', 'bead 4')
elseif b.beadsnumber==5
    legend('bead 1', 'bead 2', 'bead 3', 'bead 4', 'bead 5')
elseif b.beadsnumber==6
    legend('bead 1', 'bead 2', 'bead 3', 'bead 4', 'bead 5', 'bead 6')
end

%even frames
figure('position',[10 40 560 420]);
figure(2)
plot(b.beadsx_even,b.beadsy_even)
axis square equal
title('Even frames, beads image')
xlabel('x position (nm)')
ylabel('y position (nm)')
if b.beadsnumber==1
    legend('bead 1')
elseif b.beadsnumber==2
    legend('bead 1', 'bead 2')
elseif b.beadsnumber==3
    legend('bead 1', 'bead 2', 'bead 3')
elseif b.beadsnumber==4
    legend('bead 1', 'bead 2', 'bead 3', 'bead 4')
elseif b.beadsnumber==5
    legend('bead 1', 'bead 2', 'bead 3', 'bead 4', 'bead 5')
elseif b.beadsnumber==6
    legend('bead 1', 'bead 2', 'bead 3', 'bead 4', 'bead 5', 'bead 6')
end


%%display of individual beads and selection of OK beads for drift calculations

j_odd=0;
k_odd=0;
j_even=0;
k_even=0;
b.OK.beadsx_odd=[];
b.OK.beadsy_odd=[];
b.BAD.beadsx_odd=[];
b.BAD.beadsy_odd=[];
b.OK.beadsx_even=[];
b.OK.beadsy_even=[];
b.BAD.beadsx_even=[];
b.BAD.beadsy_even=[];

for i=1:b.beadsnumber;
    figure('position',[590 550 500 420]);
    figure(3)
    plot(b.beadsx_odd(:,i),b.beadsy_odd(:,i),'+');
    axis square equal
    title(['Odd frames, bead ', num2str(b.beadslabel(i))]);
    xlabel('x position (nm)')
    ylabel('y position (nm)')
    figure('position',[590 40 500 420]);
    figure(4)
    plot(b.beadsx_even(:,i),b.beadsy_even(:,i),'+');
    axis square equal
    title(['Even frames, bead ', num2str(b.beadslabel(i))]);
    xlabel('x position (nm)')
    ylabel('y position (nm)')
    
    figure('position',[1100 550 560 420]);
    figure(5)
    subplot(2,1,1);  plot(b.time_odd,b.beadsx_odd(:,i),'+');
    ylabel('x position (nm)')
    title(['Odd frames, bead ', num2str(b.beadslabel(i))]);
    subplot(2,1,2);  plot(b.time_odd,b.beadsy_odd(:,i),'+');
    xlabel('time (s)')
    ylabel('y position (nm)')
    
    figure('position',[1100 40 560 420]);
    figure(6)
    subplot(2,1,1);  plot(b.time_even,b.beadsx_even(:,i),'+');
    ylabel('x position (nm)')
    title(['Even frames, bead ', num2str(b.beadslabel(i))]);
    subplot(2,1,2);  plot(b.time_even,b.beadsy_even(:,i),'+');
    xlabel('time (s)')
    ylabel('y position (nm)')
    %selection of beads and labeling them as ok and bad
    flag=menu('Would you use this bead for drift calculation?','Yes','No');
    if  flag==2 %if bad bead
        k_odd=k_odd+1;
        b.BAD.beadsx_odd(:,k_odd)=b.beadsx_odd(:,i);
        b.BAD.beadsy_odd(:,k_odd)=b.beadsy_odd(:,i);
        k_even=k_even+1;
        b.BAD.indices(k_even)=b.beadslabel(i);
        b.BAD.beadsx_even(:,k_even)=b.beadsx_even(:,i);
        b.BAD.beadsy_even(:,k_even)=b.beadsy_even(:,i);
    elseif flag==1 %if OK bead
        j_odd=j_odd+1;
        b.OK.beadsx_odd(:,j_odd)=b.beadsx_odd(:,i);
        b.OK.beadsy_odd(:,j_odd)=b.beadsy_odd(:,i);
        j_even=j_even+1;
        b.OK.indices(j_even)=b.beadslabel(i);
        b.OK.beadsx_even(:,j_even)=b.beadsx_even(:,i);
        b.OK.beadsy_even(:,j_even)=b.beadsy_even(:,i);
    end
end



clear  k_odd j_odd k_even j_even i
b.OK.number=min(size(b.OK.beadsx_odd,2),size(b.OK.beadsx_even,2));
b.BAD.number=min(size(b.BAD.beadsx_odd,2),size(b.BAD.beadsx_even,2));
b.BAD.xco=[];
%
close all


% drift calculation
%calculates the center of mass of the selected beads
drift_x_odd=mean(b.OK.beadsx_odd,2);
drift_y_odd=mean(b.OK.beadsy_odd,2);
drift_x_even=mean(b.OK.beadsx_even,2);
drift_y_even=mean(b.OK.beadsy_even,2);


%
b.OK.xx_odd=b.OK.beadsx_odd;
b.OK.yy_odd=b.OK.beadsy_odd;
b.OK.xx_even=b.OK.beadsx_even;
b.OK.yy_even=b.OK.beadsy_even;
%
b.BAD.xx_odd=b.BAD.beadsx_odd;
b.BAD.yy_odd=b.BAD.beadsy_odd;
b.BAD.xx_even=b.BAD.beadsx_even;
b.BAD.yy_even=b.BAD.beadsy_even;
%


%bin number calculations
bins_number_odd=(filter_drift/0.04)+1;% filter 40 s
bins_number_even=(filter_drift/0.04)+1;% filter 40 s

x_drift_f_odd=smooth(drift_x_odd,bins_number_odd);
y_drift_f_odd=smooth(drift_y_odd,bins_number_odd);
x_drift_f_even=smooth(drift_x_even,bins_number_even);
y_drift_f_even=smooth(drift_y_even,bins_number_even);

%only time points for which the drift is calculated is used in further
%calculations
xx_drift0_odd=x_drift_f_odd(round((filter_drift/0.04)/2)+1:(length(x_drift_f_odd)-round((filter_drift/0.04)/2)));%x_drift_f_odd(bins_number_odd:bins_number_odd:length(x_drift_f_odd));
yy_drift0_odd=y_drift_f_odd(round((filter_drift/0.04)/2)+1:(length(x_drift_f_odd)-round((filter_drift/0.04)/2)));%y_drift_f_odd(bins_number_odd:bins_number_odd:length(y_drift_f_odd));
xx_drift0_even=x_drift_f_even(round((filter_drift/0.04)/2)+1:(length(x_drift_f_even)-round((filter_drift/0.04)/2)));%x_drift_f_odd(bins_number_odd:bins_number_odd:length(x_drift_f_odd));
yy_drift0_even=y_drift_f_even(round((filter_drift/0.04)/2)+1:(length(x_drift_f_even)-round((filter_drift/0.04)/2)));%y_drift_f_odd(bins_number_odd:bins_number_odd:length(y_drift_f_odd));

t_f_odd=b.time_odd(round((filter_drift/0.04)/2)+1:length(x_drift_f_odd)-round((filter_drift/0.04)/2));%(bins_number:bins_number:length(t))-t(floor(bins_number/2));
t_f_even=b.time_even(round((filter_drift/0.04)/2)+1:length(x_drift_f_even)-round((filter_drift/0.04)/2));%(bins_number:bins_number:length(t))-t(floor(bins_number/2));

Delta_driftx=mean(xx_drift0_odd-xx_drift0_even);
Delta_drifty=mean(yy_drift0_odd-yy_drift0_even);

%display of the calculated drift
figure('position',[530 550 560 420]);
figure(1)
plot(b.time_odd,drift_x_odd,'b.')
hold on
plot(t_f_odd,xx_drift0_odd,'-r')
title('Odd frames, x axis drift')
xlabel('time (s)')
ylabel('x position (nm)')
figure('position',[530 40 560 420]);
figure(2)
plot(b.time_odd,drift_y_odd,'b.')
hold on
plot(t_f_odd,yy_drift0_odd,'-r')
title('Odd frames, y axis drift')
xlabel('time (s)')
ylabel('y position (nm)')

figure('position',[1100 550 560 420]);
figure(3)
plot(b.time_even,drift_x_even,'b.')
hold on
plot(t_f_even,xx_drift0_even,'-r')
title('Even frames, x axis drift')
xlabel('time (s)')
ylabel('x position (nm)')
figure('position',[1100 40 560 420]);
figure(4)
plot(b.time_even,drift_y_even,'b.')
hold on
plot(t_f_even,yy_drift0_even,'-r')
title('Even frames, y axis drift')
xlabel('time (s)')
ylabel('y position (nm)')



%drift corrected from ok beads

for i=1:b.OK.number
    xx_drift_odd(:,i)=xx_drift0_odd;
    yy_drift_odd(:,i)=yy_drift0_odd;
    xx_drift_even(:,i)=xx_drift0_even;
    yy_drift_even(:,i)=yy_drift0_even;
end



b.OK.xco_odd=(b.OK.xx_odd(round((filter_drift/0.04)/2)+1:(length(x_drift_f_odd)-round((filter_drift/0.04)/2)),:)-xx_drift_odd);
b.OK.yco_odd=(b.OK.yy_odd(round((filter_drift/0.04)/2)+1:(length(x_drift_f_odd)-round((filter_drift/0.04)/2)),:)-yy_drift_odd);
b.OK.xco_even=(b.OK.xx_even(round((filter_drift/0.04)/2)+1:(length(x_drift_f_even)-round((filter_drift/0.04)/2)),:)-xx_drift_even);
b.OK.yco_even=(b.OK.yy_even(round((filter_drift/0.04)/2)+1:(length(x_drift_f_even)-round((filter_drift/0.04)/2)),:)-yy_drift_even);
%
mmx_odd=mean(b.OK.xco_odd,1);
mmy_odd=mean(b.OK.yco_odd,1);
mmx_even=mean(b.OK.xco_even,1);
mmy_even=mean(b.OK.yco_even,1);

for i=1:b.OK.number
    mean_x_odd(:,i) = mmx_odd(i)*ones(size(b.OK.xco_odd,1),1);
    mean_y_odd(:,i) = mmy_odd(i)*ones(size(b.OK.yco_odd,1),1);
    mean_x_even(:,i) = mmx_even(i)*ones(size(b.OK.xco_even,1),1);
    mean_y_even(:,i) = mmy_even(i)*ones(size(b.OK.yco_even,1),1);
end

b.OK.xcor_odd = b.OK.xco_odd-mean_x_odd;
b.OK.ycor_odd = b.OK.yco_odd-mean_y_odd;
b.OK.xcor_even = b.OK.xco_even-mean_x_even;
b.OK.ycor_even = b.OK.yco_even-mean_y_even;

clear xx_drift_odd xx_drift_even yy_drift_odd yy_drift_even mmx_odd mmy_odd mmx_even mmy_even mean_x_odd mean_y_odd mean_x_even mean_y_even
%drift correction of OK beads done

%drift corrected from bad beads
if b.BAD.number>=1;
    for i=1:b.BAD.number
        xx_drift_odd(:,i)=xx_drift0_odd;
        yy_drift_odd(:,i)=yy_drift0_odd;
        xx_drift_even(:,i)=xx_drift0_even;
        yy_drift_even(:,i)=yy_drift0_even;
    end
    
    b.BAD.xco_odd=(b.BAD.xx_odd(round((filter_drift/0.04)/2)+1:(length(x_drift_f_odd)-round((filter_drift/0.04)/2)),:)-xx_drift_odd);
    b.BAD.yco_odd=(b.BAD.yy_odd(round((filter_drift/0.04)/2)+1:(length(x_drift_f_odd)-round((filter_drift/0.04)/2)),:)-yy_drift_odd);
    b.BAD.xco_even=(b.BAD.xx_even(round((filter_drift/0.04)/2)+1:(length(x_drift_f_even)-round((filter_drift/0.04)/2)),:)-xx_drift_even);
    b.BAD.yco_even=(b.BAD.yy_even(round((filter_drift/0.04)/2)+1:(length(x_drift_f_even)-round((filter_drift/0.04)/2)),:)-yy_drift_even);
    
end

%above:all drift corrected. OK beads, center of amss corrected. BAD beads
%center of mass not corrected yet.



%displays the OK beads after drift correction
figure('position',[10 550 560 420]);
figure(5)
hold on
plot(b.OK.xco_odd,b.OK.yco_odd,'.')
axis square equal
title('Odd frames, beads image after drift correction')
xlabel('x position (nm)')
ylabel('y position (nm)')
if b.OK.number==1
    legend(['bead',num2str(b.OK.indices(1))])
elseif b.OK.number==2
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))])
elseif b.OK.number==3
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))])
elseif b.OK.number==4
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))],['bead',num2str(b.OK.indices(4))])
elseif b.OK.number==5
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))],['bead',num2str(b.OK.indices(4))],['bead',num2str(b.OK.indices(5))])
elseif b.OK.number==6
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))],['bead',num2str(b.OK.indices(4))],['bead',num2str(b.OK.indices(5))],['bead',num2str(b.OK.indices(6))])
end
figure('position',[10 40 560 420]);
figure(6)
hold on
plot(b.OK.xco_even,b.OK.yco_even,'.')
axis square equal
title('Even frames, beads image after drift correction')
xlabel('x position (nm)')
ylabel('y position (nm)')
if b.OK.number==1
    legend(['bead',num2str(b.OK.indices(1))])
elseif b.OK.number==2
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))])
elseif b.OK.number==3
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))])
elseif b.OK.number==4
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))],['bead',num2str(b.OK.indices(4))])
elseif b.OK.number==5
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))],['bead',num2str(b.OK.indices(4))],['bead',num2str(b.OK.indices(5))])
elseif b.OK.number==6
    legend(['bead',num2str(b.OK.indices(1))],['bead',num2str(b.OK.indices(2))],['bead',num2str(b.OK.indices(3))],['bead',num2str(b.OK.indices(4))],['bead',num2str(b.OK.indices(5))],['bead',num2str(b.OK.indices(6))])
end
pause%press any key to proceed

%all ok odd and even frames are combined
b.time(1:2:size(t_f_odd,1)+size(t_f_even,1))=t_f_odd;
b.time(2:2:size(t_f_odd,1)+size(t_f_even,1))=t_f_even;
b.OK.xcor(1:2:size(b.OK.xcor_odd,1)+size(b.OK.xcor_even,1),:)=b.OK.xcor_odd;
b.OK.xcor(2:2:size(b.OK.xcor_odd,1)+size(b.OK.xcor_even,1),:)=b.OK.xcor_even;
b.OK.ycor(1:2:size(b.OK.ycor_odd,1)+size(b.OK.ycor_even,1),:)=b.OK.ycor_odd;
b.OK.ycor(2:2:size(b.OK.ycor_odd,1)+size(b.OK.ycor_even,1),:)=b.OK.ycor_even;

%check below... I did this to take care of the drift between the odd
%and even frames. For OK beads, it is taken care of by dynamic calculation
%of the center of mass and also since we use the distance from the center
%of mass.

%all bad odd and even frames are combined
%drift between even and odd frames, corrected
if b.BAD.number>=1;
    b.BAD.xco(1:2:size(b.BAD.xco_odd,1)+size(b.BAD.xco_even,1),:)=b.BAD.xco_odd;
    b.BAD.xco(2:2:size(b.BAD.xco_odd,1)+size(b.BAD.xco_even,1),:)=b.BAD.xco_even+Delta_driftx;
    b.BAD.yco(1:2:size(b.BAD.yco_odd,1)+size(b.BAD.yco_even,1),:)=b.BAD.yco_odd;
    b.BAD.yco(2:2:size(b.BAD.yco_odd,1)+size(b.BAD.yco_even,1),:)=b.BAD.yco_even+Delta_drifty;
end


hold on
close all

% simmetry check for multiple tethers : OK beads.. the beads that fail get
% added to the BAD beads

j=0;
c.OK.indices=[];
c.REC.indices=[];
k=b.BAD.number;
for i=1:b.OK.number;
    figure('position',[10 550 560 420]);
    
    % check symmetry by calculating covariance (method by Lin Han et al.)
    lambda=eig(cov(b.OK.xcor(:,i),b.OK.ycor(:,i)));%computes the eigen value of the covariance matrix of the x and y coordinates
    lambda1=max(lambda(1),lambda(2));
    lambda2=min(lambda(1),lambda(2));
    b.OK.s(i)=sqrt(lambda1/lambda2);
    
    %
    XX=[min(b.OK.xcor(:,i)):15:max(b.OK.xcor(:,i))];
    YY=[min(b.OK.ycor(:,i)):15:max(b.OK.ycor(:,i))];
    [XHIST]=histc(b.OK.xcor(:,i),XX);
    [YHIST]=histc(b.OK.ycor(:,i),YY);
    figure('position',[530 550 560 420]);
    figure(1)
    subplot(2,1,1); bar(XX,XHIST/sum(XHIST),1);
    
    xlim([-1000 1000]);
    ylabel('probability')
    xlabel('x position (nm)')
    title(['bead',num2str(b.OK.indices(i))]);
    legend(['s', num2str(b.OK.s(i))]);
    hold off
    subplot(2,1,2); bar(YY,YHIST/sum(YHIST),1);
    
    xlim([-1000 1000]);
    ylabel('probability')
    xlabel('y position (nm)')
    hold off
    figure('position',[530 40 560 420]);
    figure(2)
    
    plot(XX,XHIST/sum(XHIST),'sb');
    hold on
    plot(YY,YHIST/sum(YHIST),'sr');
    title(['bead',num2str(b.OK.indices(i))]);
    legend(['s', num2str(b.OK.s(i))]);
    hold off
    
    % symmetry by radial histogram
    [THETA,RHO]=cart2pol(b.OK.xcor(:,i),b.OK.ycor(:,i));
    figure('position',[730 40 560 420]);
    figure(3)
    rose(THETA,60);
    hold on
    legend(['s', num2str(b.OK.s(i))]);
    hold off
    figure(4)
    hold on
    plot(b.OK.xcor(:,i),b.OK.ycor(:,i))
    axis equal
    hold off
    flag=menu('Simmetry test','Yes','No');
    if  flag==2%if bead is not symmetrical, it is added to bad beads
        k=k+1;
        b.BAD.indices(k)=b.OK.indices(i);
        b.BAD.xco(:,k)=b.OK.xcor(:,i);
        b.BAD.yco(:,k)=b.OK.ycor(:,i);
    elseif flag==1%symmetrical bead data is transfered to structure c
        j=j+1;
        c.OK.indices(j)=b.OK.indices(i);
        c.OK.beadsx(:,j)=b.OK.xcor(:,i);
        c.OK.beadsy(:,j)=b.OK.ycor(:,i);
        
    end
    clear flag XHIST YHIST XX YY
end
clear e f lambda1 lambda2
close all

b.BAD.number=size(b.BAD.xco,2);

%till here all the OK beads are good to go and all the bad beads are in BAD
%beads array. Now we select the OK area from BAD beads and reject the
%really bad beads
close all

clear j
j=0;
if b.BAD.number>=1
    for i=1:b.BAD.number;
        figure(1)
        subplot(2,1,1);  plot(b.time,b.BAD.xco(:,i),'+');
        ylabel('x position (nm)')
        title(['bead',num2str(b.BAD.indices(i))]);
        subplot(2,1,2);  plot(b.time,b.BAD.yco(:,i),'+');
        xlabel('time (s)')
        ylabel('y position (nm)')
        but= menu('Would you select some region to cut?','Yes','No');
        if but==2
        else
            h=0;
            XX=[0];
            while but == 1
                h=h+1;
                k1 = waitforbuttonpress;
                point1 = get(gca,'CurrentPoint'); % button down detected
                finalRect = rbbox;    % return figure units
                point2 = get(gca,'CurrentPoint');    % button up detected
                point1 = point1(1,1:2);              % extract x and y
                point2 = point2(1,1:2);
                p1 = min(point1,point2);             % calculate locations
                offset = abs(point1-point2);         % and dimensions
                x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
                y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
                XX=[XX, point1(1),point2(1)];
                hold on
                axis manual
                plot(x,y)
                but= menu('Any other region to cut?','Yes','No, I have finished');
            end
            XX=[XX,length(b.time)];
            clear but k1 point1 point2 p1 offset x y ccc
            ccc=[];
            for kk=1:length(XX)/2
                ccc=[ccc, find(b.time>=XX(2*kk-1) & b.time<=XX(2*kk))];
            end
            b.BAD.length(i)=length(ccc);
            for k=1:b.BAD.length(i)
                b.BAD.time(k,i)=b.time(ccc(k));
                b.BAD.xx(k,i)=b.BAD.xco(ccc(k),i);
                b.BAD.yy(k,i)=b.BAD.yco(ccc(k),i);
            end
            clear ccc
            sumx=0; sumy=0;
            for k=1:b.BAD.length(i)
                sumx=sumx+b.BAD.xx(k,i);
                sumy=sumy+b.BAD.yy(k,i);
            end
            mmx_bad(i)=sumx/b.BAD.length(i);
            mmy_bad(i)=sumy/b.BAD.length(i);
            for k=1:b.BAD.length(i)
                b.BAD.xcor(k,i)=b.BAD.xx(k,i)-mmx_bad(i);
                b.BAD.ycor(k,i)=b.BAD.yy(k,i)-mmy_bad(i);
            end
            figure(2)
            hold on
            plot(b.BAD.xcor(:,i),b.BAD.ycor(:,i))
            axis equal
            hold off
            
            %simmetry test
            % check symmetry by calculating covariance (method by Lin Han et al.)
            lambda=eig(cov(b.BAD.xcor(1:b.BAD.length(i),i),b.BAD.ycor(1:b.BAD.length(i),i)));%computes the eigen value of the covariance matrix of the x and y coordinates
            lambda1=max(lambda(1),lambda(2));
            lambda2=min(lambda(1),lambda(2));
            b.BAD.s(i)=sqrt(lambda1/lambda2);
            
            %
            XX=min(b.BAD.xcor(1:b.BAD.length(i),i)):15:max(b.BAD.xcor(1:b.BAD.length(i),i));
            YY=min(b.BAD.ycor(1:b.BAD.length(i),i)):15:max(b.BAD.ycor(1:b.BAD.length(i),i));
            [XHIST]=histc(b.BAD.xcor(1:b.BAD.length(i),i),XX);
            [YHIST]=histc(b.BAD.ycor(1:b.BAD.length(i),i),YY);
            figure('position',[530 550 560 420]);
            figure(4)
            subplot(2,1,1); bar(XX,XHIST/sum(XHIST),1);
            
            xlim([-1000 1000]);
            ylabel('probability')
            xlabel('x position (nm)')
            title(['bead',num2str(b.BAD.indices(i))]);
            legend(['s', num2str(b.BAD.s(i))]);
            hold off
            subplot(2,1,2); bar(YY,YHIST/sum(YHIST),1);
            
            xlim([-1000 1000]);
            ylabel('probability')
            xlabel('y position (nm)')
            hold off
            figure('position',[530 40 560 420]);
            figure(5)
            
            plot(XX,XHIST/sum(XHIST),'sb');
            hold on
            plot(YY,YHIST/sum(YHIST),'sr');
            title(['bead',num2str(b.BAD.indices(i))]);
            legend(['s', num2str(b.BAD.s(i))]);
            hold off
            
            % symmetry by radial histogram
            [THETA,RHO]=cart2pol(b.BAD.xcor(1:b.BAD.length(i),i),b.BAD.ycor(1:b.BAD.length(i),i));
            figure('position',[730 540 560 420]);
            figure(6)
            rose(THETA,60);
            hold on
            legend(['s', num2str(b.BAD.s(i))]);
            hold off
            
            flag=menu('Simmetry test','Yes','No');
            if  flag==2
                
            elseif flag==1
                j=j+1;
                for k=1:b.BAD.length(i)
                    c.REC.indices(j)=b.BAD.indices(i);
                    c.REC.time(k,j)=b.BAD.time(k,i);
                    c.REC.beadsx(k,j)=b.BAD.xcor(k,i);
                    c.REC.beadsy(k,j)=b.BAD.ycor(k,i);
                end
                c.REC.length(j)=b.BAD.length(i);
            end
            clear flag XHIST YHIST XX YY
            close all
        end
    end
end






%transfer the x and y to a bigger matrix so that it is same for all time
%filters

if isempty(c.OK.indices)==0
    for i=1:length(filter_rect);
        c.OK.bins(i)=bins(i); % filter times
        %bins_number_odd=(round((filter_drift/mean(diff(b.time_odd)))/2))*2+1; 
        if isempty(c.OK.indices)==0
            for z=1:size(c.OK.beadsx,2)
                for m=1:length(c.OK.beadsx)
                    c.OK.x(m,z,i)=c.OK.beadsx(m,z);
                    c.OK.y(m,z,i)=c.OK.beadsy(m,z);
                end
            end
        end
    end
    
    %find rho overall for OK beads
    for z=1:size(c.OK.x,2)
        for k=1:size(c.OK.x,1)
            c.OK.rho_overall(k,z)=sqrt((c.OK.beadsx(k,z)^2)+(c.OK.beadsy(k,z))^2);
        end
    end
    
    % rho mean and median overall
    c.OK.rho_median_overall=median(c.OK.rho_overall);%median
    c.OK.rho_mean_overall=mean(c.OK.rho_overall);%mean
    
    %make an array of average of x and y (x bar and y bar)
    for i=1:length(filter_rect)
        c.OK.bins(i)=bins(i);
        for z=1:size(c.OK.x,2)
            c.OK.x_mean(:,z,i)=smooth(c.OK.x(:,z,i),c.OK.bins(i));
            c.OK.y_mean(:,z,i)=smooth(c.OK.y(:,z,i),c.OK.bins(i));
        end
    end
    
    
    % %find rho or radius of each bead with respect to the point of attachment
    % %for that time window
    
    for i=1:length(filter_rect)
        for z=1:size(c.OK.x,2)
            for k=1:size(c.OK.x,1)
                c.OK.rho(k,z,i)=sqrt(((c.OK.x(k,z,i)-c.OK.x_mean(k,z,i))^2)+((c.OK.y(k,z,i)-c.OK.y_mean(k,z,i))^2));
                c.OK.rho_square(k,z,i)=((c.OK.x(k,z,i)-c.OK.x_mean(k,z,i))^2)+((c.OK.y(k,z,i)-c.OK.y_mean(k,z,i))^2);
            end
        end
    end
    
    %find rho average and rhosquare average in different time window
    
    for i=1:length(filter_rect)
        
        for z=1:size(c.OK.x,2)
            c.OK.rho_avg(:,z,i)=smooth(c.OK.rho(:,z,i),c.OK.bins(i));
            c.OK.rho_square_avg(:,z,i)=smooth(c.OK.rho_square(:,z,i),c.OK.bins(i));
            for k=1:size(c.OK.x,1)
                c.OK.rho_rms(k,z,i)=sqrt(c.OK.rho_square_avg(k,z,i));
            end
        end
    end
    
    %find rho median in each time window
    for i=1:length(filter_rect)
        for z=1:size(c.OK.x,2)
            for k=1:size(c.OK.x,1)
                if k>=(c.OK.bins(i)+1)/2 && k<=(size(c.OK.x,1)-(c.OK.bins(i)+1)/2)
                    c.OK.rho_median(k,z,i)=median(c.OK.rho((k+1-((c.OK.bins(i)+1)/2)):(k+((c.OK.bins(i)+1)/2)-1),z,i));
                else
                    c.OK.rho_median(k,z,i)=c.OK.rho(k,z,i);
                end
            end
        end
    end
    
    % %finds start and end of the data, sums the terms and gives the average
    clear s squaresum
    
    for k=1:size(c.OK.rho,3)
        for j=1:size(c.OK.rho,2)
            c.OK.start(1,j,1)=(round((filter_rect(length(filter_rect))/mean(diff(b.time)))/2)+1);
            c.OK.end(1,j,1)=size(c.OK.rho,1)-round((filter_rect(length(filter_rect))/mean(diff(b.time)))/2);
            
            c.OK.rho_length(1,j,k)=c.OK.end(1,j,1)-c.OK.start(1,j,1)+1;
            s=0; squaresum=0;
            for i=c.OK.start(1,j,1):c.OK.end(1,j,1)
                s=s+c.OK.rho_avg(i,j,k);
                squaresum=squaresum+c.OK.rho_rms(i,j,k);
            end
            c.OK.rho_avg_sum(1,j,k)=s;
            c.OK.rho_rms_sum(1,j,k)=squaresum;
            c.OK.rho_avg_avg(1,j,k)=c.OK.rho_avg_sum(1,j,k)/c.OK.rho_length(1,j,k);
            c.OK.rho_rms_avg(1,j,k)=c.OK.rho_rms_sum(1,j,k)/c.OK.rho_length(1,j,k);
        end
        clear s squaresum
    end
    
    % finds rho median median for the whole distribution for each time
    % window and each bead
    for k=1:size(c.OK.rho,3)
        for j=1:size(c.OK.rho,2)
            c.OK.rho_median_median(j,k)=median(c.OK.rho_median(c.OK.start(1,j,1):c.OK.end(1,j,1),j,k));
            c.OK.rho_median_mad(j,k)=mad(c.OK.rho_median(c.OK.start(1,j,1):c.OK.end(1,j,1),j,k),1);
        end
    end
    
    
    s=0;
    for k=1:size(c.OK.x,3)
        for j=1:size(c.OK.x,2)
            for i=c.OK.start(1,j,1):c.OK.end(1,j,1)
                s=s+c.OK.rho_square_avg(i,j,k);
            end
            c.OK.rho_square_avg_sum(1,j,k)=s;
            c.OK.rho_square_avg_avg(1,j,k)=c.OK.rho_square_avg_sum(1,j,k)/c.OK.rho_length(1,j,k);
            s=0;
        end
    end
    clear s t
    s=0; t=0;
    for k=1:size(c.OK.x,3)
        for j=1:size(c.OK.x,2)
            for i=c.OK.start(1,j,1):c.OK.end(1,j,1)
                s=s+(c.OK.rho_avg(i,j,k)-c.OK.rho_avg_avg(1,j,k))^2;
                t=t+(c.OK.rho_rms(i,j,k)-c.OK.rho_rms_avg(1,j,k))^2;
            end
            c.OK.rho_avg_var(1,j,k)=s/(c.OK.rho_length(1,j,k)-1);
            c.OK.rho_avg_std(1,j,k)=sqrt(s/(c.OK.rho_length(1,j,k)-1));
            c.OK.rho_rms_var(1,j,k)=t/(c.OK.rho_length(1,j,k)-1);
            c.OK.rho_rms_std(1,j,k)=sqrt(t/(c.OK.rho_length(1,j,k)-1));
            c.OK.rho_stdt(1,j,k)=sqrt(c.OK.rho_square_avg_avg(1,j,k)-(c.OK.rho_rms_avg(1,j,k)^2));
            s=0; t=0;
        end
    end
    clear s t
    
    %storing the time windows in the structure c
    c.REC.filter_rect=filter_rect;
    c.OK.filter_rect=filter_rect;
    
    %plots for OK beads
    %plot 1
    clear x y figura
    for j=1:size(c.OK.x,2)
        for k=1:length(filter_rect)
            x(k)=c.OK.rho_avg_avg(1,j,k);
            y(k)=c.OK.rho_avg_std(1,j,k);
        end
        figura=figure(1);
        subplot(2,2,1);
        plot(filter_rect,x,'--r*');
        hold on
        subplot(2,2,2);
        loglog(filter_rect,x,'--g^');
        hold on
        subplot(2,2,3);
        plot(filter_rect,y,'--rs');
        hold on
        subplot(2,2,4);
        loglog(filter_rect,y,'--g+');
        hold on
        if j==size(c.OK.x,2)
            hold off;
            saveas(figura,[pathname,'figures\','bead ',' controlOK1'],'fig')
            close
        end
    end
    clear x y figura
    
    
    %plot 2
    
    for j=1:size(c.OK.x,2)
        for k=1:length(filter_rect)
            x(k)=c.OK.rho_rms_avg(1,j,k);
            y(k)=c.OK.rho_rms_std(1,j,k);
        end
        figura=figure(2);
        subplot(2,2,1);
        plot(filter_rect,x,'--r*');
        hold on
        subplot(2,2,2);
        loglog(filter_rect,x,'--g^');
        hold on
        subplot(2,2,3);
        plot(filter_rect,y,'--rs');
        hold on
        subplot(2,2,4);
        loglog(filter_rect,y,'--g+');
        hold on
        if j==size(c.OK.x,2)
            hold off;
            saveas(figura,[pathname,'figures\','bead ',' controlOK2'],'fig')
            close
        end
    end
    clear x y figura
    
    
    %plot 3
    
    for j=1:size(c.OK.x,2)
        for k=1:length(filter_rect)
            x(k)=c.OK.rho_stdt(1,j,k);
        end
        figura=figure(3);
        subplot(2,1,1);
        plot(filter_rect,x,'--r*');
        hold on
        subplot(2,1,2);
        loglog(filter_rect,x,'--g^');
        hold on
        if j==size(c.OK.x,2)
            hold off;
            saveas(figura,[pathname,'figures\','bead ',' controlOK3'],'fig')
            close
        end
    end
    clear x y figura
end



%transfer the x and y to a bigger matrix so that it is same for all time
%filters
if isempty(c.REC.indices)==0
    for i=1:length(filter_rect);
        c.REC.bins(i)=bins(i); % filter times
        
        for z=1:size(c.REC.beadsx,2)
            for m=1:length(c.REC.beadsx)
                c.REC.x(m,z,i)=c.REC.beadsx(m,z);
                c.REC.y(m,z,i)=c.REC.beadsy(m,z);
            end
            
        end
    end
    
    %REC beads
    %find mean x and mean y overall
    for j=1:size(c.REC.x,2)
        c.REC.meanx_overall(j)=mean(c.REC.beadsx(1:c.REC.length(j),j));
        c.REC.meany_overall(j)=mean(c.REC.beadsy(1:c.REC.length(j),j));
    end
    
    %find rho overall for REC beads
    
    for z=1:size(c.REC.x,2)
        for k=1:size(c.REC.x,1)
            c.REC.rho_overall(k,z)=sqrt((c.REC.beadsx(k,z)-c.REC.meanx_overall(z))^2+(c.REC.beadsy(k,z)-c.REC.meany_overall(j))^2);
        end
    end
    
    % rho mean and median overall
    for i=1:size(c.REC.x,2)
        c.REC.rho_median_overall(i)=median(c.REC.rho_overall(1:c.REC.length(i),i));%median
        c.REC.rho_mean_overall(i)=mean(c.REC.rho_overall(1:c.REC.length(i),i));%mean
    end
    %make an array of average of x and y (x bar and y bar)
    
    for i=1:length(filter_rect)
        c.REC.bins(i)=bins(i);
        for z=1:size(c.REC.x,2)
            c.REC.x_mean(:,z,i)=smooth(c.REC.x(:,z,i),c.REC.bins(i));
            c.REC.y_mean(:,z,i)=smooth(c.REC.y(:,z,i),c.REC.bins(i));
        end
    end
    
    %
    % %find rho or radius of each bead with respect to the point of attachment
    % %for that time window
    
    for i=1:length(filter_rect)
        
        for z=1:size(c.REC.x,2)
            for k=1:size(c.REC.x,1)
                c.REC.rho(k,z,i)=sqrt(((c.REC.x(k,z,i)-c.REC.x_mean(k,z,i))^2)+((c.REC.y(k,z,i)-c.REC.y_mean(k,z,i))^2));
                c.REC.rho_square(k,z,i)=((c.REC.x(k,z,i)-c.REC.x_mean(k,z,i))^2)+((c.REC.y(k,z,i)-c.REC.y_mean(k,z,i))^2);
            end
        end
    end
    
    %find rho average and rhosquare average in different time window
    
    for i=1:length(filter_rect)
       
        for z=1:size(c.REC.x,2)
            c.REC.rho_avg(:,z,i)=smooth(c.REC.rho(:,z,i),c.REC.bins(i));
            c.REC.rho_square_avg(:,z,i)=smooth(c.REC.rho_square(:,z,i),c.REC.bins(i));
            for k=1:size(c.REC.x,1)
                c.REC.rho_rms(k,z,i)=sqrt(c.REC.rho_square_avg(k,z,i));
            end
        end
    end
    
    %find rho median in each time window
    for i=1:length(filter_rect)
        for z=1:size(c.REC.x,2)
            for k=1:size(c.REC.x,1)
                if k>=(c.REC.bins(i)+1)/2 && k<=(size(c.REC.x,1)-(c.REC.bins(i)+1)/2)
                    c.REC.rho_median(k,z,i)=median(c.REC.rho((k+1-((c.REC.bins(i)+1)/2)):(k+((c.REC.bins(i)+1)/2)-1),z,i));
                else
                    c.REC.rho_median(k,z,i)=c.REC.rho(k,z,i);
                end
            end
        end
    end
    
    
    % %finds start and end of the data, sums the terms and gives the average
    clear s squaresum
    
    for k=1:size(c.REC.rho,3)
        for j=1:size(c.REC.rho,2)
            c.REC.start(1,j,1)=(round((filter_rect(length(filter_rect))/mean(diff(b.time)))/2)+1);
            c.REC.end(1,j,1)=c.REC.length(j)-round((filter_rect(length(filter_rect))/mean(diff(b.time)))/2);
            
            c.REC.rho_length(1,j,k)=c.REC.end(1,j,1)-c.REC.start(1,j,1)+1;
            s=0; squaresum=0;
            for i=c.REC.start(1,j,1):c.REC.end(1,j,1)
                s=s+c.REC.rho_avg(i,j,k);
                squaresum=squaresum+c.REC.rho_rms(i,j,k);
            end
            c.REC.rho_avg_sum(1,j,k)=s;
            c.REC.rho_rms_sum(1,j,k)=squaresum;
            c.REC.rho_avg_avg(1,j,k)=c.REC.rho_avg_sum(1,j,k)/c.REC.rho_length(1,j,k);
            c.REC.rho_rms_avg(1,j,k)=c.REC.rho_rms_sum(1,j,k)/c.REC.rho_length(1,j,k);
        end
        clear s squaresum
    end
    
    %finds rho median median
    for k=1:size(c.REC.rho,3)
        for j=1:size(c.REC.rho,2)
            c.REC.rho_median_median(j,k)=median(c.REC.rho_median((c.REC.start(1,j,1):c.REC.end(1,j,1)),j,k));
            c.REC.rho_median_mad(j,k)=mad(c.REC.rho_median(c.REC.start(1,j,1):c.REC.end(1,j,1),j,k),1);
        end
    end
    
    s=0;
    for k=1:size(c.REC.x,3)
        for j=1:size(c.REC.x,2)
            for i=c.REC.start(1,j,1):c.REC.end(1,j,1)
                s=s+c.REC.rho_square_avg(i,j,k);
            end
            c.REC.rho_square_avg_sum(1,j,k)=s;
            c.REC.rho_square_avg_avg(1,j,k)=c.REC.rho_square_avg_sum(1,j,k)/c.REC.rho_length(1,j,k);
            s=0;
        end
    end
    clear s t
    s=0; t=0;
    for k=1:size(c.REC.x,3)
        for j=1:size(c.REC.x,2)
            for i=c.REC.start(1,j,1):c.REC.end(1,j,1)
                s=s+(c.REC.rho_avg(i,j,k)-c.REC.rho_avg_avg(1,j,k))^2;
                t=t+(c.REC.rho_rms(i,j,k)-c.REC.rho_rms_avg(1,j,k))^2;
                
            end
            c.REC.rho_avg_var(1,j,k)=s/(c.REC.rho_length(1,j,k)-1);
            c.REC.rho_avg_std(1,j,k)=sqrt(s/(c.REC.rho_length(1,j,k)-1));
            c.REC.rho_rms_var(1,j,k)=t/(c.REC.rho_length(1,j,k)-1);
            c.REC.rho_rms_std(1,j,k)=sqrt(t/(c.REC.rho_length(1,j,k)-1));
            c.REC.rho_stdt(1,j,k)=sqrt(c.REC.rho_square_avg_avg(1,j,k)-(c.REC.rho_rms_avg(1,j,k)^2));
            s=0; t=0;
        end
    end
    clear s t
    
    
    
    %plots for REC beads
    %plot 1
    clear x y
    for j=1:size(c.REC.x,2)
        for k=1:length(filter_rect)
            x(k)=c.REC.rho_avg_avg(1,j,k);
            y(k)=c.REC.rho_avg_std(1,j,k);
        end
        figura=figure(1);
        subplot(2,2,1);
        plot(filter_rect,x,'--r*');
        hold on
        subplot(2,2,2);
        loglog(filter_rect,x,'--g^');
        hold on
        subplot(2,2,3);
        plot(filter_rect,y,'--rs');
        hold on
        subplot(2,2,4);
        loglog(filter_rect,y,'--g+');
        hold on
        if j==size(c.REC.x,2)
            hold off;
            saveas(figura,[pathname,'figures\','bead ',' controlREC1'],'fig')
            close
        end
    end
    clear x y figura
    
    
    %plot 2
    
    for j=1:size(c.REC.x,2)
        for k=1:length(filter_rect)
            x(k)=c.REC.rho_rms_avg(1,j,k);
            y(k)=c.REC.rho_rms_std(1,j,k);
        end
        figura=figure(2);
        subplot(2,2,1);
        plot(filter_rect,x,'--r*');
        hold on
        subplot(2,2,2);
        loglog(filter_rect,x,'--g^');
        hold on
        subplot(2,2,3);
        plot(filter_rect,y,'--rs');
        hold on
        subplot(2,2,4);
        loglog(filter_rect,y,'--g+');
        hold on
        if j==size(c.REC.x,2)
            hold off;
            saveas(figura,[pathname,'figures\','bead ',' controlREC2'],'fig')
            close
        end
    end
    clear x y figura
    
    
    %plot 3
    
    for j=1:size(c.REC.x,2)
        for k=1:length(filter_rect)
            x(k)=c.REC.rho_stdt(1,j,k);
        end
        figura=figure(3);
        subplot(2,1,1);
        plot(filter_rect,x,'--r*');
        hold on
        subplot(2,1,2);
        loglog(filter_rect,x,'--g^');
        hold on
        if j==size(c.REC.x,2)
            hold off;
            saveas(figura,[pathname,'figures\','bead ',' controlREC3'],'fig')
            close
        end
    end
    clear x y figura
    
end

%storage of data in structure c
c.SOURCE=[pathname,filename(1:length(filename)-4)];%path of starting text file
save([pathname,filename(1:length(filename)-4),'.mat'], 'c')%structure with all the analysed data stored
close all
clear all
clc
done='done'