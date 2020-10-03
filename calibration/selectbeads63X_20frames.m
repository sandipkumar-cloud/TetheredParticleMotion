%SELECTION OF BEAD AND STORAGE OF COORDINATES, TIME AND LENGTH OF THE TRACE

%TPM data analysis CONTROL $interlaced files$ $50Hz$ 
%correction of the pixel size
%incorporation of better symmetry test.. covariance matrix (lin  han et al.)
%http://en.wikipedia.org/wiki/Covariance_matrix
%above link explains why we take the ratio of the square root of the eigen
%values(square root of the eigen values represent the lengths of the major and minor axes)
%plot radial histogram

clear all

close all
% filters variable
filter_drift=20/50;% 20 frames
filter_drift1=10;% If there are no stuck beads, the 10s window average is subtracted from the data to remove drift



%%%%%%%%%%%%
% chose file to analyze
ext='.txt';
[filename, pathname] = uigetfile({'*.txt';'*.*'},'File Selector');

%remove figures folder and make a new one
if exist([pathname,'figures'],'dir')
    rmdir([pathname,'figures'],'s');
end

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
legend('show')

%even frames
figure('position',[10 40 560 420]);
figure(2)
plot(b.beadsx_even,b.beadsy_even)
axis square equal
title('Even frames, beads image')
xlabel('x position (nm)')
ylabel('y position (nm)')
legend('show')


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

if isempty(b.OK.beadsx_odd)
    bins_number_odd1=(filter_drift1*25)+1;
    bins_number_even1=(filter_drift1*25)+1;
    t_f_odd=b.time_odd(round((filter_drift1*25)/2)+1:(size(b.BAD.beadsx_odd,1)-round((filter_drift1*25)/2)));
    t_f_even=b.time_even(round((filter_drift1*25)/2)+1:(size(b.BAD.beadsx_odd,1)-round((filter_drift1*25)/2)));
    b.time(1:2:size(t_f_odd,1)+size(t_f_even,1))=t_f_odd;
    b.time(2:2:size(t_f_odd,1)+size(t_f_even,1))=t_f_even;
    for i=1:b.BAD.number
        for j=1:size(b.BAD.beadsx_odd,1)
            xx_odd(j)=b.BAD.beadsx_odd(j,i);
            yy_odd(j)=b.BAD.beadsy_odd(j,i);
            xx_even(j)=b.BAD.beadsx_even(j,i);
            yy_even(j)=b.BAD.beadsy_even(j,i);
        end
        
        drift_x_odd=smooth(xx_odd,bins_number_odd1);
        drift_y_odd=smooth(yy_odd,bins_number_odd1);
        drift_x_even=smooth(xx_even,bins_number_even1);
        drift_y_even=smooth(yy_even,bins_number_even1);
        
        badxodd=xx_odd.'-drift_x_odd;
        badyodd=yy_odd.'-drift_y_odd;
        badxeven=xx_even.'-drift_x_even;
        badyeven=yy_even.'-drift_y_even;
        k=1;
       for j=round((filter_drift1*25)/2)+1:(size(b.BAD.beadsx_odd,1)-round((filter_drift1*25)/2))
        b.BAD.xco_odd(k,i)=badxodd(j,1);
        b.BAD.yco_odd(k,i)=badyodd(j,1);
        b.BAD.xco_even(k,i)=badxeven(j,1);
        b.BAD.yco_even(k,i)=badyeven(j,1);
        k=k+1;
       end
    end
else


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
bins_number_odd=(filter_drift*25)+1;
bins_number_even=(filter_drift*25)+1;


x_drift_f_odd=smooth(drift_x_odd,bins_number_odd);
y_drift_f_odd=smooth(drift_y_odd,bins_number_odd);
x_drift_f_even=smooth(drift_x_even,bins_number_even);
y_drift_f_even=smooth(drift_y_even,bins_number_even);

%only time points for which the drift is calculated is used in further
%calculations
xx_drift0_odd=x_drift_f_odd(round((filter_drift*25)/2)+1:(length(x_drift_f_odd)-round((filter_drift*25)/2)));
yy_drift0_odd=y_drift_f_odd(round((filter_drift*25)/2)+1:(length(x_drift_f_odd)-round((filter_drift*25)/2)));
xx_drift0_even=x_drift_f_even(round((filter_drift*25)/2)+1:(length(x_drift_f_even)-round((filter_drift*25)/2)));
yy_drift0_even=y_drift_f_even(round((filter_drift*25)/2)+1:(length(x_drift_f_even)-round((filter_drift*25)/2)));

t_f_odd=b.time_odd(round((filter_drift*25)/2)+1:length(x_drift_f_odd)-round((filter_drift*25)/2));
t_f_even=b.time_even(round((filter_drift*25)/2)+1:length(x_drift_f_even)-round((filter_drift*25)/2));



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



b.OK.xco_odd=(b.OK.xx_odd(round((filter_drift*25)/2)+1:(length(x_drift_f_odd)-round((filter_drift*25)/2)),:)-xx_drift_odd);
b.OK.yco_odd=(b.OK.yy_odd(round((filter_drift*25)/2)+1:(length(x_drift_f_odd)-round((filter_drift*25)/2)),:)-yy_drift_odd);
b.OK.xco_even=(b.OK.xx_even(round((filter_drift*25)/2)+1:(length(x_drift_f_even)-round((filter_drift*25)/2)),:)-xx_drift_even);
b.OK.yco_even=(b.OK.yy_even(round((filter_drift*25)/2)+1:(length(x_drift_f_even)-round((filter_drift*25)/2)),:)-yy_drift_even);
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
    
    b.BAD.xco_odd=(b.BAD.xx_odd(round((filter_drift*25)/2)+1:(length(x_drift_f_odd)-round((filter_drift*25)/2)),:)-xx_drift_odd);
    b.BAD.yco_odd=(b.BAD.yy_odd(round((filter_drift*25)/2)+1:(length(x_drift_f_odd)-round((filter_drift*25)/2)),:)-yy_drift_odd);
    b.BAD.xco_even=(b.BAD.xx_even(round((filter_drift*25)/2)+1:(length(x_drift_f_even)-round((filter_drift*25)/2)),:)-xx_drift_even);
    b.BAD.yco_even=(b.BAD.yy_even(round((filter_drift*25)/2)+1:(length(x_drift_f_even)-round((filter_drift*25)/2)),:)-yy_drift_even);
    
end
end
%above:all drift corrected. OK beads, center of mass corrected. BAD beads
%center of mass not corrected yet.


if isempty(b.OK.beadsx_odd)
else
%displays the OK beads after drift correction
figure('position',[10 550 560 420]);
figure(5)
hold on
plot(b.OK.xco_odd,b.OK.yco_odd,'.')
axis square equal
title('Odd frames, beads image after drift correction')
xlabel('x position (nm)')
ylabel('y position (nm)')
legend('show')
figure('position',[10 40 560 420]);
figure(6)
hold on
plot(b.OK.xco_even,b.OK.yco_even,'.')
axis square equal
title('Even frames, beads image after drift correction')
xlabel('x position (nm)')
ylabel('y position (nm)')
legend('show')
pause%press any key to proceed

%all ok odd and even frames are combined
b.time(1:2:size(t_f_odd,1)+size(t_f_even,1))=t_f_odd;
b.time(2:2:size(t_f_odd,1)+size(t_f_even,1))=t_f_even;
b.OK.xcor(1:2:size(b.OK.xcor_odd,1)+size(b.OK.xcor_even,1),:)=b.OK.xcor_odd;
b.OK.xcor(2:2:size(b.OK.xcor_odd,1)+size(b.OK.xcor_even,1),:)=b.OK.xcor_even;
b.OK.ycor(1:2:size(b.OK.ycor_odd,1)+size(b.OK.ycor_even,1),:)=b.OK.ycor_odd;
b.OK.ycor(2:2:size(b.OK.ycor_odd,1)+size(b.OK.ycor_even,1),:)=b.OK.ycor_even;

end

%check below... I did this to take care of the drift between the odd
%and even frames. For OK beads, it is taken care of by dynamic calculation
%of the center of mass and also since we use the distance from the center
%of mass.

%all bad odd and even frames are combined
%drift between even and odd frames, corrected
if b.BAD.number>=1;
    b.BAD.xco(1:2:size(b.BAD.xco_odd,1)+size(b.BAD.xco_even,1),:)=b.BAD.xco_odd;
    b.BAD.xco(2:2:size(b.BAD.xco_odd,1)+size(b.BAD.xco_even,1),:)=b.BAD.xco_even;
    b.BAD.yco(1:2:size(b.BAD.yco_odd,1)+size(b.BAD.yco_even,1),:)=b.BAD.yco_odd;
    b.BAD.yco(2:2:size(b.BAD.yco_odd,1)+size(b.BAD.yco_even,1),:)=b.BAD.yco_even;
end


hold on
close all

% simmetry check for multiple tethers : OK beads.. the beads that fail get
% added to the BAD beads

j=0;
c.OK.indices=[];
c.REC.indices=[];
k=b.BAD.number;
if isempty(b.OK.beadsx_odd)
else
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
    figure('position',[730 540 560 420]);
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
        c.OK.time=b.time;
        c.OK.indices(j)=b.OK.indices(i);
        c.OK.beadsx(:,j)=b.OK.xcor(:,i);
        c.OK.beadsy(:,j)=b.OK.ycor(:,i);
        c.OK.length(j)=size(c.OK.beadsx,1);
    end
    close all
    clear flag XHIST YHIST XX YY
end
clear e f lambda1 lambda2
close all
end

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
        close all
    end
end
%symmetry test done
%% 
%storage of data in structure c
c.SOURCE=[pathname,filename(1:length(filename)-4)];%path of starting text file
save([pathname,filename(1:length(filename)-4),'.mat'], 'c')%structure with all the analysed data stored
close all
clear all
clc
done='done'